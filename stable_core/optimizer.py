"""
RFOptimizer: Physical Stability-Constrained RF Parameter Optimizer.

Replaces traditional flat-potential (xi = 1.0) solving by enforcing physical
stability constraints (PTBL, Coupled D-Q, LCBI) via penalty functions.
Maximizes bunch lengthening / Touschek lifetime while strictly guaranteeing STABLE beam.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Callable, Dict, Optional, Tuple, Any, Sequence

import numpy as np
from scipy.optimize import minimize, differential_evolution

from stable_core.types import (
    MachineParameters,
    StabilityStatus,
    InstabilityResult,
)
from stable_core.classifier import InstabilityClassifier

logger = logging.getLogger(__name__)


@dataclass
class OptimizationResult:
    """Result of stability-constrained RF optimization."""
    success: bool
    optimal_params: Dict[str, float]
    bunch_length_factor: float
    xi: float
    status: StabilityStatus
    growth_rate: float
    nyquist_margin: float
    iterations: int
    message: str
    instability_result: Optional[InstabilityResult] = None
    convergence_history: list[dict[str, Any]] = field(default_factory=list)


class RFOptimizer:
    """
    Physical Stability-Constrained RF Parameter Optimizer.

    Traditional flat-potential solving targets xi = 1.0 (zero first and second
    potential derivatives at bunch center), where the synchrotron frequency
    omega_s0 vanishes, making the beam highly susceptible to PTBL or FMCI.

    This optimizer treats beam dynamic stability as an essential hard penalty constraint:
    Objective = - (Bunch Lengthening Factor) + Penalty(StabilityStatus, GrowthRate, NyquistMargin)
    """

    def __init__(
        self,
        machine: Optional[MachineParameters] = None,
        classifier: Optional[InstabilityClassifier] = None,
        penalty_unstable: float = 1e4,
        penalty_margin_weight: float = 50.0,
    ) -> None:
        self.machine = machine or MachineParameters()
        self.classifier = classifier or InstabilityClassifier(self.machine)
        self.penalty_unstable = penalty_unstable
        self.penalty_margin_weight = penalty_margin_weight

    def optimize_harmonic_cavity_detuning(
        self,
        I0_mA: float,
        initial_detuning_khz: Optional[float] = None,
        detuning_bounds_khz: Tuple[float, float] = (20.0, 150.0),
        method: str = "COBYLA",
        evaluate_func: Optional[Callable[[float, float], InstabilityResult]] = None,
        target_margin: float = 0.20,
    ) -> OptimizationResult:
        """
        Find optimal harmonic cavity detuning Delta f_hc maximizing bunch lengthening
        while guaranteeing physical stability.

        Args:
            I0_mA: Beam current in mA.
            initial_detuning_khz: Starting guess for Delta f_hc in kHz. If None, derived from xi ~ 0.95.
            detuning_bounds_khz: (min_khz, max_khz) search interval.
            method: Optimization algorithm ('COBYLA', 'Nelder-Mead', 'Powell', or 'diff_evolution').
            evaluate_func: Optional custom evaluator (I0_mA, df_khz) -> InstabilityResult.
                           Defaults to self.classifier.evaluate_semianalytical.
            target_margin: Desired minimum Nyquist / stability distance margin.

        Returns:
            OptimizationResult containing optimal detuning, bunch lengthening factor, and stability status.
        """
        if evaluate_func is None:
            evaluate_func = lambda i0, df: self.classifier.evaluate_semianalytical(i0, df)

        # Initial guess from analytical formula if not given
        if initial_detuning_khz is None:
            # Flat potential detuning:
            # xi = I0 * R_hc / (n_hc * V_mc * (R/Q)_hc * ...) ~ approx
            # Search detuning corresponding to xi ~ 0.95
            nom_df = self.machine.get_nominal_hc_detuning(I0_mA)
            initial_detuning_khz = float(np.clip(nom_df, detuning_bounds_khz[0], detuning_bounds_khz[1]))

        history: list[dict[str, Any]] = []

        def objective(x: Sequence[float]) -> float:
            df_khz = float(x[0])

            # Bound check penalty
            if df_khz < detuning_bounds_khz[0] or df_khz > detuning_bounds_khz[1]:
                return self.penalty_unstable * 2.0

            res = evaluate_func(I0_mA, df_khz)
            xi = self.classifier.compute_flat_potential_parameter(I0_mA, df_khz)

            # Base cost: negative bunch length factor (we want to maximize bunch length)
            cost = -res.bunch_length_factor

            # Hard penalty if beam is physically unstable
            if not res.is_stable or res.growth_rate > 0.0:
                cost += self.penalty_unstable + 1000.0 * max(0.0, res.growth_rate)
                cost += 500.0 * (int(res.status.value) + 1)
            else:
                # Soft penalty for low stability margin
                if res.nyquist_margin < target_margin:
                    cost += self.penalty_margin_weight * (target_margin - res.nyquist_margin)

            history.append({
                "df_khz": df_khz,
                "xi": xi,
                "bunch_length_factor": res.bunch_length_factor,
                "status": res.status.name,
                "growth_rate": res.growth_rate,
                "nyquist_margin": res.nyquist_margin,
                "cost": cost,
            })
            return cost

        # Run optimization
        if method.lower() == "diff_evolution":
            res_opt = differential_evolution(
                objective,
                bounds=[detuning_bounds_khz],
                seed=42,
                maxiter=50,
                tol=1e-3,
            )
            opt_df = float(res_opt.x[0])
            success = res_opt.success
            message = str(res_opt.message)
            n_iter = res_opt.nit
        else:
            opt_bounds = [detuning_bounds_khz] if method.upper() in ["COBYLA", "POWELL", "L-BFGS-B"] else None
            # COBYLA uses constraints for bounds
            constraints = []
            if method.upper() == "COBYLA":
                constraints = [
                    {"type": "ineq", "fun": lambda x: x[0] - detuning_bounds_khz[0]},
                    {"type": "ineq", "fun": lambda x: detuning_bounds_khz[1] - x[0]},
                ]

            res_opt = minimize(
                objective,
                x0=np.array([initial_detuning_khz]),
                method=method,
                bounds=opt_bounds if method.upper() != "COBYLA" else None,
                constraints=constraints if method.upper() == "COBYLA" else (),
                options={"maxiter": 100, "disp": False},
            )
            opt_df = float(res_opt.x[0])
            success = res_opt.success
            message = str(res_opt.message)
            n_iter = res_opt.nit if hasattr(res_opt, "nit") else res_opt.nfev

        # Final verification at optimal point
        final_eval = evaluate_func(I0_mA, opt_df)
        final_xi = self.classifier.compute_flat_potential_parameter(I0_mA, opt_df)

        return OptimizationResult(
            success=success and final_eval.is_stable,
            optimal_params={"detuning_khz": opt_df, "I0_mA": I0_mA},
            bunch_length_factor=final_eval.bunch_length_factor,
            xi=final_xi,
            status=final_eval.status,
            growth_rate=final_eval.growth_rate,
            nyquist_margin=final_eval.nyquist_margin,
            iterations=int(n_iter),
            message=message,
            instability_result=final_eval,
            convergence_history=history,
        )

    def optimize_multivariable_rf(
        self,
        I0_mA: float,
        variable_bounds: Dict[str, Tuple[float, float]],
        initial_guess: Optional[Dict[str, float]] = None,
        method: str = "COBYLA",
    ) -> OptimizationResult:
        """
        Multivariable RF parameter optimization (e.g., HC detuning and Direct Feedback Gain).

        Variables supported:
            - 'detuning_khz': HC detuning in kHz
            - 'K_p': Direct RF feedback gain
        """
        keys = list(variable_bounds.keys())
        bounds = [variable_bounds[k] for k in keys]

        if initial_guess is None:
            x0 = [0.5 * (b[0] + b[1]) for b in bounds]
        else:
            x0 = [initial_guess.get(k, 0.5 * (variable_bounds[k][0] + variable_bounds[k][1])) for k in keys]

        history: list[dict[str, Any]] = []

        def objective(x: Sequence[float]) -> float:
            p_dict = {k: float(v) for k, v in zip(keys, x)}

            # Check bounds
            for k, (b_lo, b_hi) in variable_bounds.items():
                if p_dict[k] < b_lo or p_dict[k] > b_hi:
                    return self.penalty_unstable * 2.0

            df_khz = p_dict.get("detuning_khz", 40.0)
            kp = p_dict.get("K_p", 1.0)

            # Evaluate instability
            res = self.classifier.evaluate_semianalytical(I0_mA, df_khz, K_p=kp)
            xi = self.classifier.compute_flat_potential_parameter(I0_mA, df_khz)

            cost = -res.bunch_length_factor
            if not res.is_stable or res.growth_rate > 0.0:
                cost += self.penalty_unstable + 1000.0 * max(0.0, res.growth_rate)
            else:
                if res.nyquist_margin < 0.20:
                    cost += self.penalty_margin_weight * (0.20 - res.nyquist_margin)

            history.append({**p_dict, "cost": cost, "status": res.status.name, "growth_rate": res.growth_rate})
            return cost

        constraints = []
        if method.upper() == "COBYLA":
            for i, (b_lo, b_hi) in enumerate(bounds):
                constraints.append({"type": "ineq", "fun": lambda x, i=i, b_lo=b_lo: x[i] - b_lo})
                constraints.append({"type": "ineq", "fun": lambda x, i=i, b_hi=b_hi: b_hi - x[i]})

        res_opt = minimize(
            objective,
            x0=np.array(x0),
            method=method,
            bounds=bounds if method.upper() != "COBYLA" else None,
            constraints=constraints if method.upper() == "COBYLA" else (),
            options={"maxiter": 150, "disp": False},
        )

        opt_p = {k: float(v) for k, v in zip(keys, res_opt.x)}
        df_khz = opt_p.get("detuning_khz", 40.0)
        kp = opt_p.get("K_p", 1.0)
        final_eval = self.classifier.evaluate_semianalytical(I0_mA, df_khz, K_p=kp)
        final_xi = self.classifier.compute_flat_potential_parameter(I0_mA, df_khz)

        return OptimizationResult(
            success=res_opt.success and final_eval.is_stable,
            optimal_params=opt_p,
            bunch_length_factor=final_eval.bunch_length_factor,
            xi=final_xi,
            status=final_eval.status,
            growth_rate=final_eval.growth_rate,
            nyquist_margin=final_eval.nyquist_margin,
            iterations=int(res_opt.nit if hasattr(res_opt, "nit") else res_opt.nfev),
            message=str(res_opt.message),
            instability_result=final_eval,
            convergence_history=history,
        )
