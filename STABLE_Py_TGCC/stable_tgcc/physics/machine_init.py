"""
Exact Python port of MATLAB STABLE initialization routines:
  - calc_bunch_form_factor_eq.m
  - detune_HC_calc.m
  - computeOptimalMCDetune (STABLE_App.m)
  - alignHOMToWorstCaseSideband (STABLE_App.m)
  - generatePattern (STABLE_App.m)
  - TruncatedGaussian.m
  - machine.m
"""

from __future__ import annotations
from typing import Dict, Any, Tuple
import math
import cmath
import numpy as np
from scipy.special import erf, erfinv


def truncated_gaussian(sigma: float, limits: Tuple[float, float], shape: Tuple[int, ...]) -> np.ndarray:
    """Exact port of TruncatedGaussian.m using inverse error function."""
    a, b = float(limits[0]), float(limits[1])
    cdf_a = 0.5 * (1.0 + erf(a / (sigma * math.sqrt(2.0))))
    cdf_b = 0.5 * (1.0 + erf(b / (sigma * math.sqrt(2.0))))
    u = np.random.uniform(cdf_a, cdf_b, size=shape)
    return (sigma * math.sqrt(2.0) * erfinv(2.0 * u - 1.0)).astype(np.float32)


def generate_fill_pattern(h: int, mode: str = "Uniform (100%)", custom_expr: str = "") -> np.ndarray:
    """Generate bucket fill pattern (1=filled, 0=empty) matching STABLE_App.generatePattern."""
    h = int(h)
    m = mode.strip()
    if m == "80% filling":
        pattern = np.ones(h, dtype=np.int32)
        n_groups = h // 40
        for ii in range(1, n_groups + 1):
            start = 40 * ii - 8
            end = min(h, 40 * ii)
            pattern[start:end] = 0
        return pattern
    elif m == "90% filling":
        pattern = np.ones(h, dtype=np.int32)
        n_groups = h // 40
        for ii in range(1, n_groups + 1):
            start = 40 * ii - 4
            end = min(h, 40 * ii)
            pattern[start:end] = 0
        return pattern
    elif m == "Single bunch":
        pattern = np.zeros(h, dtype=np.int32)
        pattern[0] = 1
        return pattern
    elif m == "Custom" and custom_expr.strip():
        try:
            arr = np.asarray(eval(custom_expr, {"np": np, "ones": np.ones, "zeros": np.zeros, "h": h}), dtype=np.int32).ravel()
            if len(arr) == h:
                return (arr != 0).astype(np.int32)
        except Exception:
            pass
    return np.ones(h, dtype=np.int32)


def calc_bunch_form_factor_eq(
    mode: str,
    I0: float,
    V_mc: float,
    U0: float,
    E0: float,
    C: float,
    h: int,
    alpha_c: float,
    sigma_e0: float,
    n_hc: int,
    R_hc: float,
    Q_hc: float,
    param_val: float,
    use_form_factor: bool = True,
) -> Dict[str, Any]:
    """
    Self-consistent Fokker-Planck potential-well equilibrium solver for passive harmonic cavity
    bunch form factor F_n(sigma_tau) & flat-potential detuning (exact port of calc_bunch_form_factor_eq.m).
    """
    cspeed = 299792458.0
    T0 = C / cspeed
    f0 = 1.0 / T0
    f_rf = h * f0
    w_rf = 2.0 * math.pi * f_rf
    f_hc = n_hc * f_rf
    wre_hc = 2.0 * math.pi * f_hc

    k_fp_sq = 1.0 / (n_hc ** 2) - (1.0 / max(n_hc ** 2 - 1.0, 1.0)) * (U0 / max(V_mc, 1.0)) ** 2
    k_fp = math.sqrt(max(k_fp_sq, 1e-8))

    Qs0 = math.sqrt(max(0.0, alpha_c / (2.0 * math.pi) * h * math.sqrt(max(0.0, (V_mc / E0) ** 2 - (U0 / E0) ** 2))))
    sigma_t0 = alpha_c * sigma_e0 / max(2.0 * math.pi * Qs0, 1e-12) * T0

    out: Dict[str, Any] = {"sigma_t0": sigma_t0, "k_fp": k_fp}

    if mode.lower() == "detune_for_xi":
        xi_target = float(param_val)
        if I0 <= 1e-9 or R_hc <= 0 or Q_hc <= 0 or k_fp <= 0:
            out.update({
                "fre_shift_hc": 0.0,
                "psi_hc_deg": 0.0,
                "F_hc": 1.0,
                "F_mc": 1.0,
                "sigma_t_eq": sigma_t0,
                "xi_eff": 0.0,
            })
            return out

        if not use_form_factor:
            cos_psi = (xi_target * k_fp * V_mc) / (2.0 * I0 * R_hc)
            if cos_psi >= 1.0:
                psi = 0.0
                df_hc = 0.0
            else:
                psi = math.acos(max(-1.0, cos_psi))
                df_hc = math.tan(psi) * f_hc / (2.0 * Q_hc)
            out.update({
                "fre_shift_hc": df_hc,
                "psi_hc_deg": psi * 180.0 / math.pi,
                "F_hc": 1.0,
                "F_mc": 1.0,
                "sigma_t_eq": sigma_t0,
                "xi_eff": (2.0 * I0 * R_hc * math.cos(psi)) / (k_fp * V_mc),
            })
            return out

        tau_max = min(0.45 / f_rf, max(12.0 * sigma_t0, 320e-12))
        tau = np.linspace(-tau_max, tau_max, 401)
        dtau = float(tau[1] - tau[0])

        Fn = 0.92
        F1 = 0.995
        th1 = 0.0
        thn = 0.0
        psi = 0.0
        sig_t = sigma_t0

        for it in range(1, 26):
            cos_psi = (xi_target * k_fp * V_mc) / max(2.0 * Fn * I0 * R_hc, 1e-9)
            psi = 0.0 if cos_psi >= 0.999999 else math.acos(max(0.0, cos_psi))

            arg_s = (U0 + 2.0 * (Fn ** 2) * I0 * R_hc * (math.cos(psi) ** 2)) / max(F1 * V_mc, 1.0)
            arg_s = max(-0.999, min(0.999, arg_s))
            phis = math.pi - math.asin(arg_s) - th1

            Fn_c = Fn * cmath.exp(-1j * thn)
            E_gain = (
                V_mc * np.sin(phis + w_rf * tau)
                - 2.0 * I0 * R_hc * math.cos(psi) * np.real(cmath.exp(1j * psi) * Fn_c * np.exp(1j * n_hc * w_rf * tau))
                - U0
            )
            U_pot = np.concatenate(([0.0], -np.cumsum(0.5 * (E_gain[:-1] + E_gain[1:])) * dtau)) / (T0 * E0)
            U_pot = U_pot - np.min(U_pot)

            exponent = np.clip(U_pot / max(alpha_c * (sigma_e0 ** 2), 1e-18), 0.0, 100.0)
            rho = np.exp(-exponent)
            rho = rho / max(float(np.sum(rho)) * dtau, 1e-30)

            t_mean = float(np.sum(tau * rho)) * dtau
            tau_c = tau - t_mean
            sig_t = math.sqrt(float(np.sum((tau_c ** 2) * rho)) * dtau)

            F1_new = complex(np.sum(rho * np.exp(-1j * w_rf * tau_c)) * dtau)
            Fn_new = complex(np.sum(rho * np.exp(-1j * n_hc * w_rf * tau_c)) * dtau)

            err = abs(abs(Fn_new) - Fn)
            Fn = 0.5 * abs(Fn_new) + 0.5 * Fn
            F1 = 0.5 * abs(F1_new) + 0.5 * F1
            th1 = -cmath.phase(F1_new)
            thn = -cmath.phase(Fn_new)

            if err < 1e-6 and it >= 4:
                break

        df_hc = math.tan(psi) * f_hc / (2.0 * Q_hc)
        out.update({
            "fre_shift_hc": df_hc,
            "psi_hc_deg": psi * 180.0 / math.pi,
            "F_hc": Fn,
            "F_mc": F1,
            "sigma_t_eq": sig_t,
            "xi_eff": (2.0 * Fn * I0 * R_hc * math.cos(psi)) / (k_fp * V_mc),
        })
        return out

    # Mode = 'eval'
    fre_shift_hc = float(param_val)
    wr_hc = wre_hc + 2.0 * math.pi * fre_shift_hc
    det_angle_hc = math.pi - math.atan(Q_hc * (wr_hc / wre_hc - wre_hc / wr_hc))
    psi_hc = math.pi - det_angle_hc

    if (not use_form_factor) or I0 <= 1e-9 or R_hc <= 0 or Q_hc <= 0:
        arg_s = max(-0.9999, min(0.9999, (U0 + 2.0 * I0 * R_hc * (math.cos(det_angle_hc) ** 2)) / max(V_mc, 1.0)))
        out.update({
            "F_hc": 1.0,
            "F_mc": 1.0,
            "theta_hc": 0.0,
            "theta_mc": 0.0,
            "sigma_t_eq": sigma_t0,
            "det_angle_hc": det_angle_hc,
            "psi_hc_deg": psi_hc * 180.0 / math.pi,
            "U_hc_real": 2.0 * I0 * R_hc * (math.cos(psi_hc) ** 2),
            "fais_mc_whc": math.pi - math.asin(arg_s),
            "V_hc_load_0": -2.0 * I0 * R_hc * math.cos(det_angle_hc) * cmath.exp(1j * (math.pi - det_angle_hc)),
            "xi_eff": (2.0 * I0 * R_hc * abs(math.cos(psi_hc))) / max(k_fp * V_mc, 1.0),
        })
        return out

    tau_max = min(0.45 / f_rf, max(12.0 * sigma_t0, 320e-12))
    tau = np.linspace(-tau_max, tau_max, 401)
    dtau = float(tau[1] - tau[0])

    Fn_c: complex = 0.90 + 0j
    F1_c: complex = 0.994 + 0j
    sig_t = sigma_t0
    converged = False
    phis_cent = math.pi - math.asin(min(0.999, U0 / max(V_mc, 1.0)))

    for it in range(1, 36):
        arg_s = (U0 + 2.0 * (abs(Fn_c) ** 2) * I0 * R_hc * (math.cos(psi_hc) ** 2)) / max(abs(F1_c) * V_mc, 1.0)
        arg_s = max(-0.999, min(0.999, arg_s))
        phis = math.pi - math.asin(arg_s) + cmath.phase(F1_c)

        E_gain = (
            V_mc * np.sin(phis + w_rf * tau)
            - 2.0 * I0 * R_hc * math.cos(psi_hc) * np.real(cmath.exp(1j * psi_hc) * Fn_c * np.exp(1j * n_hc * w_rf * tau))
            - U0
        )
        U_pot = np.concatenate(([0.0], -np.cumsum(0.5 * (E_gain[:-1] + E_gain[1:])) * dtau)) / (T0 * E0)
        U_pot = U_pot - np.min(U_pot)

        exponent = np.clip(U_pot / max(alpha_c * (sigma_e0 ** 2), 1e-18), 0.0, 100.0)
        rho = np.exp(-exponent)
        rho = rho / max(float(np.sum(rho)) * dtau, 1e-30)

        t_mean = float(np.sum(tau * rho)) * dtau
        tau_c = tau - t_mean
        sig_t = math.sqrt(float(np.sum((tau_c ** 2) * rho)) * dtau)

        F1_new = complex(np.sum(rho * np.exp(-1j * w_rf * tau_c)) * dtau)
        Fn_new = complex(np.sum(rho * np.exp(-1j * n_hc * w_rf * tau_c)) * dtau)

        err = abs(Fn_new - Fn_c)
        Fn_c = 0.35 * Fn_new + 0.65 * Fn_c
        F1_c = 0.35 * F1_new + 0.65 * F1_c
        phis_cent = phis + w_rf * t_mean

        if err < 1e-6 and it >= 5:
            converged = True
            break

    xi_raw = (2.0 * abs(Fn_c) * I0 * R_hc * abs(math.cos(psi_hc))) / max(k_fp * V_mc, 1.0)
    if (not converged) or phis_cent > (math.pi - 0.03) or phis_cent < (0.5 * math.pi) or (xi_raw > 1.06 and sig_t < 1.8 * sigma_t0):
        fp_ref = calc_bunch_form_factor_eq(
            "detune_for_xi", I0, V_mc, U0, E0, C, h, alpha_c, sigma_e0, n_hc, R_hc, Q_hc, 1.04, True
        )
        Fn_mag = float(fp_ref["F_hc"])
        F1_mag = float(fp_ref["F_mc"])
        th_hc = 0.0
        sig_t = float(fp_ref["sigma_t_eq"])
        arg_s = max(-0.999, min(0.999, (U0 + 2.0 * (Fn_mag ** 2) * I0 * R_hc * (math.cos(psi_hc) ** 2)) / max(F1_mag * V_mc, 1.0)))
        phis_cent = math.pi - math.asin(arg_s)
    else:
        Fn_mag = abs(Fn_c)
        F1_mag = abs(F1_c)
        th_hc = -cmath.phase(Fn_c)

    out.update({
        "F_hc": Fn_mag,
        "F_mc": F1_mag,
        "theta_hc": th_hc,
        "sigma_t_eq": sig_t,
        "det_angle_hc": det_angle_hc,
        "psi_hc_deg": psi_hc * 180.0 / math.pi,
        "U_hc_real": 2.0 * (Fn_mag ** 2) * I0 * R_hc * (math.cos(psi_hc) ** 2),
        "fais_mc_whc": phis_cent,
        "V_hc_load_0": -2.0 * Fn_mag * I0 * R_hc * math.cos(det_angle_hc) * cmath.exp(1j * (math.pi - det_angle_hc - th_hc)),
        "xi_eff": (2.0 * Fn_mag * I0 * R_hc * abs(math.cos(psi_hc))) / max(k_fp * V_mc, 1.0),
    })
    return out


def detune_hc_calc(
    I0: float,
    n_hc: int,
    C: float,
    h: int,
    U0: float,
    V_mc: float,
    R_hc: float,
    Q_hc: float,
    E0: float = 2.75e9,
    alpha_c: float = 1.07e-4,
    sigma_e0: float = 1.06e-3,
    use_form_factor: bool = True,
    xi_target: float = 1.0,
) -> float:
    """Exact port of detune_HC_calc.m."""
    res = calc_bunch_form_factor_eq(
        "detune_for_xi", I0, V_mc, U0, E0, C, h, alpha_c, sigma_e0, n_hc, R_hc, Q_hc, xi_target, use_form_factor
    )
    return float(res["fre_shift_hc"])


def compute_optimal_mc_detune(
    C: float,
    h: int,
    I0: float,
    V_mc: float,
    U0: float,
    R_mc: float,
    Q_mc: float,
    R_hc: float,
    Q_hc: float,
    fre_shift_hc: float,
    n_hc: int,
) -> float:
    """
    Exact Robinson reactive beam-loading compensation detuning (port of STABLE_App.computeOptimalMCDetune):
      tan(psi_mc,opt) = -2 * I0 * R_mc * sin(theta_s) / V_mc
    """
    cspeed = 299792458.0
    f_rf = h / (C / cspeed)
    f_hc = n_hc * f_rf
    psi_hc = math.atan(2.0 * Q_hc * fre_shift_hc / max(f_hc, 1.0))
    delta_hc = 2.0 * I0 * R_hc * (math.cos(psi_hc) ** 2)
    ratio = min(0.9999, max(-0.9999, (U0 + delta_hc) / max(V_mc, 1.0)))
    theta_s = math.acos(ratio)
    tan_psi_mc = -(2.0 * I0 * R_mc * math.sin(theta_s)) / max(V_mc, 1.0)
    return tan_psi_mc * f_rf / (2.0 * max(Q_mc, 1.0))


def align_hom_to_worst_sideband(
    fre_hom_Hz: np.ndarray,
    C: float,
    h: int,
    E0: float,
    U0: float,
    V_mc: float,
    alpha_c: float,
) -> np.ndarray:
    """
    Align each HOM frequency to its nearest positive coupled-bunch upper synchrotron sideband
    f_hom_aligned = (p * h + mu + Q_s) * f0 (port of STABLE_App.alignHOMToWorstCaseSideband).
    """
    cspeed = 299792458.0
    f0 = cspeed / C
    Qs = math.sqrt(max(0.0, alpha_c / (2.0 * math.pi) * h * math.sqrt(max(0.0, (V_mc / E0) ** 2 - (U0 / E0) ** 2))))
    aligned = np.array(fre_hom_Hz, dtype=np.float64, copy=True)
    for idx in range(len(aligned)):
        k_rev = round(aligned[idx] / f0 - Qs)
        aligned[idx] = (k_rev + Qs) * f0
    return aligned


def build_half_machine(
    C: float,
    I0: float,
    U0: float,
    E0: float,
    tau_s: float,
    tau_z: float,
    sigma_t0: float,
    sigma_e0: float,
    alpha_c: float,
    h: int,
    V_mc: float,
    n_hc: int,
    R_hc: float,
    Q_hc: float,
    fillrate: float,
    fre_shift_hc: float,
    Q_mc: float,
    R_mc: float,
    fre_shift_mc: float,
    use_form_factor: bool = False,
) -> Dict[str, Any]:
    """
    Exact Python port of machine.m.
    Returns dictionary `HALF` containing all derived physical quantities and initial cavity phasors.
    """
    cspeed = 299792458.0
    HALF: Dict[str, Any] = {}
    HALF["C"] = float(C)
    HALF["R"] = float(C) / (2.0 * math.pi)
    HALF["I0"] = float(I0)
    HALF["h"] = int(h)
    HALF["U0"] = float(U0)
    HALF["E0"] = float(E0)
    HALF["sigma_t0"] = float(sigma_t0)
    HALF["sigma_e0"] = float(sigma_e0)
    HALF["alpha_c"] = float(alpha_c)

    # Fokker-Planck fluctuation-dissipation consistency (worst-case longer damping time)
    if abs(tau_s - tau_z) > 1e-9:
        tau_long = max(float(tau_s), float(tau_z))
        tau_s = tau_long
        tau_z = tau_long

    HALF["tau_s"] = float(tau_s)
    HALF["tau_z"] = float(tau_z)
    HALF["V_mc"] = float(V_mc)
    HALF["n_hc"] = int(n_hc)
    HALF["R_hc"] = float(R_hc)
    HALF["Q_hc"] = float(Q_hc)
    HALF["fillrate"] = float(fillrate)
    HALF["fre_shift_hc"] = float(fre_shift_hc)
    HALF["R_mc"] = float(R_mc)
    HALF["Q_mc"] = float(Q_mc)
    HALF["fre_shift_mc"] = float(fre_shift_mc)
    HALF["use_form_factor"] = bool(use_form_factor)

    gamma = E0 / 0.51099906 / 1e6
    beta = math.sqrt(1.0 - 1.0 / (gamma ** 2))
    HALF["gamma"] = gamma
    HALF["beta"] = beta

    T0 = HALF["C"] / cspeed / beta
    Tb = T0 / h
    f_rf = h / T0
    w_rf = f_rf * 2.0 * math.pi
    HALF["T0"] = T0
    HALF["Tb"] = Tb
    HALF["f_rf"] = f_rf
    HALF["w_rf"] = w_rf

    fre_hc = f_rf * n_hc
    wre_hc = fre_hc * 2.0 * math.pi
    qc = T0 * I0 / h / max(fillrate, 1e-9)
    wr_hc = w_rf * n_hc + fre_shift_hc * 2.0 * math.pi
    angle_hc = wr_hc * Tb
    HALF["fre_hc"] = fre_hc
    HALF["wre_hc"] = wre_hc
    HALF["qc"] = qc
    HALF["wr_hc"] = wr_hc
    HALF["angle_hc"] = angle_hc

    fre_mc = f_rf
    wre_mc = fre_mc * 2.0 * math.pi
    wr_mc = wre_mc + fre_shift_mc * 2.0 * math.pi
    angle_mc = wr_mc * Tb
    HALF["fre_mc"] = fre_mc
    HALF["wre_mc"] = wre_mc
    HALF["wr_mc"] = wr_mc
    HALF["angle_mc"] = angle_mc

    HALF["lambda_rads"] = T0 / tau_s
    HALF["lambda_radz"] = T0 / tau_z
    HALF["fais_nat"] = math.pi - math.asin(min(0.9999, U0 / max(V_mc, 1.0)))

    det_angle_hc = math.pi - math.atan(Q_hc * (wr_hc / wre_hc - wre_hc / wr_hc))
    HALF["det_angle_hc"] = det_angle_hc
    HALF["psi_hc_deg"] = (math.pi - det_angle_hc) * 180.0 / math.pi

    if use_form_factor:
        fp_eq = calc_bunch_form_factor_eq(
            "eval", I0, V_mc, U0, E0, C, h, alpha_c, sigma_e0, n_hc, R_hc, Q_hc, fre_shift_hc, True
        )
        HALF["F_hc_eq"] = float(fp_eq["F_hc"])
        HALF["F_mc_eq"] = float(fp_eq["F_mc"])
        HALF["sigma_t_eq"] = float(fp_eq["sigma_t_eq"])
        HALF["xi_eff"] = float(fp_eq["xi_eff"])
    else:
        fp_eq = calc_bunch_form_factor_eq(
            "eval", I0, V_mc, U0, E0, C, h, alpha_c, sigma_e0, n_hc, R_hc, Q_hc, fre_shift_hc, False
        )
        HALF["F_hc_eq"] = 1.0
        HALF["F_mc_eq"] = 1.0
        HALF["sigma_t_eq"] = sigma_t0
        HALF["xi_eff"] = float(fp_eq["xi_eff"])

    # Point-bunch turn-0 initial cavity phasor & MC reference phase (avoids turn-0 dipole shock)
    arg_mc = min(0.9999, max(-0.9999, (U0 + 2.0 * I0 * R_hc * (math.cos(det_angle_hc) ** 2)) / max(V_mc, 1.0)))
    HALF["fais_mc_whc"] = math.pi - math.asin(arg_mc)
    HALF["V_hc_load_0"] = -2.0 * I0 * R_hc * math.cos(det_angle_hc) * cmath.exp(1j * (math.pi - det_angle_hc))

    det_angle_mc = math.atan(Q_mc * (wr_mc / wre_mc - wre_mc / wr_mc))
    HALF["det_angle_mc"] = det_angle_mc
    HALF["V_mc_load_0"] = 2.0 * I0 * R_mc * math.cos(det_angle_mc) * cmath.exp(1j * det_angle_mc)

    HALF["Vrf_ideal"] = -V_mc * cmath.exp(-1j * (math.pi / 2.0 - HALF["fais_mc_whc"]))
    HALF["Vg_mc_init"] = HALF["Vrf_ideal"] - HALF["V_mc_load_0"]
    HALF["Vg_mc_track_0"] = HALF["Vg_mc_init"]
    HALF["Ig_mc_0"] = HALF["Vg_mc_track_0"] / (2.0 * R_mc * math.cos(det_angle_mc)) * cmath.exp(-1j * det_angle_mc)
    HALF["Vg_mc_0"] = 2.0 * math.pi * (R_mc / Q_mc) * HALF["Ig_mc_0"]

    HALF["drift_coef"] = alpha_c * T0 * (sigma_e0 / sigma_t0)
    HALF["kick_coef"] = 1.0 / (E0 * sigma_e0)
    HALF["rfcoef1"] = HALF["kick_coef"] * V_mc
    HALF["rfcoef2"] = w_rf * sigma_t0
    HALF["radampcoef"] = 2.0 * HALF["lambda_rads"]
    HALF["quanexcoef"] = 2.0 * math.sqrt(HALF["lambda_radz"])
    HALF["ploss"] = HALF["kick_coef"] * U0

    HALF["VbImagFactor_hc"] = 1.0 / math.sqrt(4.0 * (Q_hc ** 2) - 1.0)
    HALF["rot_coef_hc"] = math.sqrt(1.0 - 1.0 / ((2.0 * Q_hc) ** 2))
    HALF["VbImagFactor_mc"] = 1.0 / math.sqrt(4.0 * (Q_mc ** 2) - 1.0)
    HALF["rot_coef_mc"] = math.sqrt(1.0 - 1.0 / ((2.0 * Q_mc) ** 2))

    return HALF
