"""
SafetyBoundaryExtractor: Automatic Safety Boundary Extractor and Equation Fitter.

Extracts critical boundary points between STABLE and UNSTABLE regions in 2D/3D
parameter spaces (e.g. R1/Q1 vs R2/Q2, detuning vs I0) via bisection or grid contouring,
and fits closed-form engineering safety boundary equations:
    e.g., detuning_khz > a * I0_mA + b, or R1/Q1 + alpha * R2/Q2 < C.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Callable, List, Optional, Tuple, Sequence

import numpy as np

from stable_core.types import StabilityStatus, InstabilityResult
from stable_core.cascade_scanner import ScanGridResult

logger = logging.getLogger(__name__)


@dataclass
class BoundaryFitResult:
    """Fitted safety boundary model and statistics."""
    x_name: str
    y_name: str
    equation_str: str
    inequality_str: str
    coefficients: np.ndarray       # highest degree first: [a_n, ..., a_1, a_0]
    degree: int
    r2_score: float
    rmse: float
    critical_x: np.ndarray
    critical_y: np.ndarray
    safe_region: str               # 'above' (y > f(x)) or 'below' (y < f(x))

    def evaluate_boundary(self, x: np.ndarray | float) -> np.ndarray | float:
        """Calculate the boundary y-value for given x."""
        poly = np.poly1d(self.coefficients)
        return poly(x)

    def is_safe(self, x: float, y: float) -> bool:
        """Check if an operating point (x, y) satisfies the safe boundary equation."""
        y_bound = float(self.evaluate_boundary(x))
        if self.safe_region == "above":
            return y > y_bound
        else:
            return y < y_bound


class SafetyBoundaryExtractor:
    """
    Extracts critical transition coordinates and fits engineering safety boundary formulas.
    """

    def __init__(self, target_status: StabilityStatus = StabilityStatus.STABLE) -> None:
        self.target_status = target_status

    def bisection_find_threshold(
        self,
        eval_func: Callable[[float, float], bool],
        x_val: float,
        y_min: float,
        y_max: float,
        tol: float = 1e-3,
        max_iter: int = 30,
    ) -> Optional[float]:
        """
        1D adaptive bisection search for the critical transition coordinate y*(x_val)
        where eval_func(x_val, y) flips between True (STABLE) and False (UNSTABLE).
        """
        s_min = eval_func(x_val, y_min)
        s_max = eval_func(x_val, y_max)

        if s_min == s_max:
            # No transition detected within interval
            return None

        lo, hi = y_min, y_max
        for _ in range(max_iter):
            mid = 0.5 * (lo + hi)
            if (hi - lo) < tol:
                return mid

            s_mid = eval_func(x_val, mid)
            if s_mid == s_min:
                lo = mid
            else:
                hi = mid

        return 0.5 * (lo + hi)

    def extract_from_evaluator(
        self,
        eval_func: Callable[[float, float], bool],
        x_name: str,
        x_values: Sequence[float],
        y_name: str,
        y_bounds: Tuple[float, float],
        degree: int = 1,
        safe_region: str = "above",
    ) -> BoundaryFitResult:
        """
        Extract safety boundary by running bisection search across an array of x_values.

        Args:
            eval_func: Callable(x, y) -> bool returning True if STABLE, False if UNSTABLE.
            x_name: Label for x-axis parameter (e.g. 'I0_mA' or 'R1_over_Q1').
            x_values: Grid points along x to search for transition.
            y_name: Label for y-axis parameter (e.g. 'detuning_khz' or 'R2_over_Q2').
            y_bounds: (y_min, y_max) search range.
            degree: Polynomial degree (1 for linear, 2 for parabolic).
            safe_region: 'above' (safe if y > f(x)) or 'below' (safe if y < f(x)).

        Returns:
            BoundaryFitResult with fitted parameters and equation.
        """
        crit_x = []
        crit_y = []

        for x in x_values:
            y_star = self.bisection_find_threshold(
                eval_func=eval_func,
                x_val=x,
                y_min=y_bounds[0],
                y_max=y_bounds[1],
            )
            if y_star is not None:
                crit_x.append(x)
                crit_y.append(y_star)

        if len(crit_x) < degree + 1:
            raise ValueError(f"Insufficient critical points ({len(crit_x)}) found to fit polynomial of degree {degree}.")

        cx = np.array(crit_x)
        cy = np.array(crit_y)

        return self._fit_polynomial(cx, cy, x_name, y_name, degree=degree, safe_region=safe_region)

    def extract_from_grid(
        self,
        grid_result: ScanGridResult,
        safe_region: str = "above",
        degree: int = 1,
    ) -> BoundaryFitResult:
        """
        Extract critical boundary from a 2D ScanGridResult by finding contour transitions.
        """
        stable_mask = grid_result.get_stable_mask()
        xs = grid_result.x_values
        ys = grid_result.y_values

        crit_x = []
        crit_y = []

        # Scan along columns (each x coordinate) to find transition in y
        for i, x in enumerate(xs):
            col = stable_mask[:, i]
            # Detect where col flips from 0 to 1 or 1 to 0
            diff = np.diff(col.astype(int))
            transitions = np.where(diff != 0)[0]
            for idx in transitions:
                # Linear interpolation between ys[idx] and ys[idx+1]
                y_trans = 0.5 * (ys[idx] + ys[idx + 1])
                crit_x.append(x)
                crit_y.append(y_trans)

        if len(crit_x) < degree + 1:
            raise ValueError(f"Found only {len(crit_x)} boundary points on grid; cannot fit degree {degree}.")

        return self._fit_polynomial(
            np.array(crit_x),
            np.array(crit_y),
            x_name=grid_result.x_name,
            y_name=grid_result.y_name,
            degree=degree,
            safe_region=safe_region,
        )

    def _fit_polynomial(
        self,
        cx: np.ndarray,
        cy: np.ndarray,
        x_name: str,
        y_name: str,
        degree: int,
        safe_region: str,
    ) -> BoundaryFitResult:
        """Fit polynomial and construct equation strings."""
        coeffs = np.polyfit(cx, cy, deg=degree)
        poly = np.poly1d(coeffs)
        y_pred = poly(cx)

        # Statistics
        ss_res = np.sum((cy - y_pred) ** 2)
        ss_tot = np.sum((cy - np.mean(cy)) ** 2)
        r2 = 1.0 - (ss_res / (ss_tot + 1e-12))
        rmse = float(np.sqrt(np.mean((cy - y_pred) ** 2)))

        # Format equation
        op = ">" if safe_region == "above" else "<"
        if degree == 1:
            slope, intercept = coeffs[0], coeffs[1]
            sign = "+" if intercept >= 0 else "-"
            eq_str = f"{y_name} = {slope:.4e} * {x_name} {sign} {abs(intercept):.4e}"
            ineq_str = f"{y_name} {op} {slope:.4e} * {x_name} {sign} {abs(intercept):.4e}"
        elif degree == 2:
            a, b, c = coeffs[0], coeffs[1], coeffs[2]
            eq_str = f"{y_name} = {a:.4e} * {x_name}^2 + {b:.4e} * {x_name} + {c:.4e}"
            ineq_str = f"{y_name} {op} {a:.4e} * {x_name}^2 + {b:.4e} * {x_name} + {c:.4e}"
        else:
            terms = [f"{coeffs[i]:.4e} * {x_name}^{degree - i}" for i in range(degree)]
            terms.append(f"{coeffs[-1]:.4e}")
            eq_str = f"{y_name} = " + " + ".join(terms)
            ineq_str = f"{y_name} {op} " + " + ".join(terms)

        return BoundaryFitResult(
            x_name=x_name,
            y_name=y_name,
            equation_str=eq_str,
            inequality_str=ineq_str,
            coefficients=coeffs,
            degree=degree,
            r2_score=float(r2),
            rmse=rmse,
            critical_x=cx,
            critical_y=cy,
            safe_region=safe_region,
        )
