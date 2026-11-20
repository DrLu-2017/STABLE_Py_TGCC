"""Instability evaluation pipeline and physical classifier.

Combines:
1. LLRF Nyquist stability / phase-gain margin screening
2. Spatial Discrete Fourier Transform (Spatial DFT) for multi-bunch mode decomposition
3. Hilbert transform envelope extraction and exponential growth rate fitting
4. Physical multi-mode state categorization (STABLE, PTBL, COUPLED_DQ, LCBI, LOST)
"""

from __future__ import annotations

import numpy as np
from scipy.signal import hilbert
from typing import Tuple, Optional, Dict, Any

from stable_core.types import MachineParameters, StabilityStatus, InstabilityResult


def _cav_impedance_and_deriv(w: complex, R: float, Q: float, wr: float) -> Tuple[complex, complex]:
    """Analytic resonator impedance Z(w) and derivative dZ/dw (Eq. 2 of T. He et al., PRAB 2024)."""
    if R <= 0.0 or Q <= 0.0 or wr <= 0.0:
        return 0.0 + 0.0j, 0.0 + 0.0j
    u = Q * (wr / w - w / wr)
    den = 1.0 + 1j * u
    if abs(den) < 1e-14:
        den = 1e-14 + 0.0j
    Z = R / den
    du_dw = -Q * (wr / (w * w) + 1.0 / wr)
    dZ_dw = -R * 1j * du_dw / (den * den)
    return Z, dZ_dw


def _eval_coupled_robinson_dispersion(
    Om: complex,
    ws2: float,
    tau_z_s: float,
    xi: float,
    n_harmonic: int,
    R1_eff: float,
    Q1_eff: float,
    wr1: float,
    w10: float,
    Rn: float,
    Qn: float,
    wrn: float,
    wn0: float,
) -> Tuple[complex, complex]:
    """Evaluate exact coupled dipole Robinson dispersion g(Omega) and derivative g'(Omega) (Eq. 8)."""
    if R1_eff > 0.0 and Q1_eff > 0.0 and wr1 > 0.0:
        Z1p, dZ1p = _cav_impedance_and_deriv(w10 + Om, R1_eff, Q1_eff, wr1)
        Z1m, dZ1m = _cav_impedance_and_deriv(-w10 + Om, R1_eff, Q1_eff, wr1)
        dZ1 = Z1p - Z1m
        ddZ1 = dZ1p - dZ1m
    else:
        dZ1 = ddZ1 = 0.0 + 0.0j

    Znp, dZnp = _cav_impedance_and_deriv(wn0 + Om, Rn, Qn, wrn)
    Znm, dZnm = _cav_impedance_and_deriv(-wn0 + Om, Rn, Qn, wrn)
    dZn = Znp - Znm
    ddZn = dZnp - dZnm

    g = (
        Om * Om
        + 1j * (2.0 / tau_z_s) * Om
        - ws2
        - 1j * xi * (dZ1 + n_harmonic * dZn)
    )
    dg = (
        2.0 * Om
        + 1j * (2.0 / tau_z_s)
        - 1j * xi * (ddZ1 + n_harmonic * ddZn)
    )
    return g, dg


def _solve_dmode_at_detune(
    df_rn_Hz: float,
    I0: float,
    ws0_sq: float,
    xi: float,
    tau_z_s: float,
    w0: float,
    h: int,
    R1_eff: float,
    Q1_eff: float,
    df_r1_Hz: float,
    n_harmonic: int,
    Rn: float,
    Qn: float,
) -> Dict[str, float]:
    """Solve He et al. (PRAB 2024) closed-form D-mode formula + exact Newton refinement."""
    w10 = h * w0
    wn0 = n_harmonic * h * w0
    dw_r1 = 2.0 * np.pi * df_r1_Hz
    dw_rn = 2.0 * np.pi * df_rn_Hz
    wr1 = w10 + dw_r1
    wrn = wn0 + dw_rn

    if abs(dw_rn) < 2.0 * np.pi * 1.0 or I0 <= 1e-9 or Rn <= 0.0 or Qn <= 0.0:
        return {
            "Omega_r": 0.0,
            "Omega_i": -1.0 / tau_z_s,
            "delta_w1": 0.0,
            "Omega_r_he": 0.0,
            "Omega_i_he": -1.0 / tau_z_s,
            "Omega_C_r": np.sqrt(max(ws0_sq, 0.0)),
            "Omega_C_i": -1.0 / tau_z_s,
        }

    # Equilibrium HHC synchrotron frequency shift (Eq. 5 of T. He et al., PRAB 2024)
    u_n0 = Qn * (wrn / wn0 - wn0 / wrn)
    Im_Zn0 = -Rn * u_n0 / (1.0 + u_n0 * u_n0)
    ws2 = ws0_sq + 2.0 * xi * n_harmonic * Im_Zn0
    Im_Z10 = 0.0
    if R1_eff > 0.0 and Q1_eff > 0.0:
        u_10 = Q1_eff * (wr1 / w10 - w10 / wr1)
        Im_Z10 = -R1_eff * u_10 / (1.0 + u_10 * u_10)
    ws2_C = max(ws2 + 2.0 * xi * Im_Z10, 1.0)

    # 1. Main-cavity phase shifts & zeta1, b1, k1 (Eqs. 8-14)
    if R1_eff > 0.0 and Q1_eff > 0.0:
        psi_minus = np.arctan(2.0 * Q1_eff * (dw_r1 - abs(dw_rn)) / wr1)
        psi_plus = np.arctan(2.0 * Q1_eff * (dw_r1 + abs(dw_rn)) / wr1)
        zeta1 = -0.5 * xi * R1_eff * (np.sin(2.0 * psi_minus) + np.sin(2.0 * psi_plus))
        b1 = 0.5 * xi * R1_eff * (np.cos(2.0 * psi_minus) - np.cos(2.0 * psi_plus))
        k1 = -(2.0 * xi * R1_eff * Q1_eff / wr1) * (
            (np.cos(psi_minus) ** 2) * np.cos(2.0 * psi_minus)
            - (np.cos(psi_plus) ** 2) * np.cos(2.0 * psi_plus)
        )
    else:
        zeta1 = b1 = k1 = 0.0

    # 2. Quadratic equation for delta_w1 & closed-form Omega_r, Omega_i (Eqs. 10-14 / Appendix B)
    dw_abs = abs(dw_rn)
    disc = 1.0
    if dw_rn >= 0.0:
        B1 = (
            (dw_rn / 4.0)
            + ((zeta1 - ws2) / (4.0 * dw_rn))
            - (xi * n_harmonic * Rn * wrn) / (16.0 * Qn * (dw_rn ** 2))
        )
        B2 = (xi * n_harmonic * wrn * Rn) / (4.0 * dw_rn * Qn)
        disc = B1 * B1 - B2
        if disc >= 0.0 and B1 > 0.0:
            delta_w1_he = B1 - np.sqrt(disc)
            if abs(delta_w1_he) < 1e-6:
                delta_w1_he = 1e-6
            Omega_r_he = dw_rn - delta_w1_he
            bn = (xi * n_harmonic * Rn * (wrn ** 2)) / (4.0 * (Qn ** 2) * (delta_w1_he ** 2))
            kn = (xi * n_harmonic * Rn * wrn) / (2.0 * Qn * (delta_w1_he ** 2))
            denom_he = 2.0 * Omega_r_he - kn - k1
            if abs(denom_he) < 1e-9:
                denom_he = -1e-9 if denom_he <= 0.0 else 1e-9
            Omega_i_he = (bn + b1 - (2.0 * Omega_r_he / tau_z_s)) / denom_he
            Om_init_he = complex(max(Omega_r_he, 10.0), np.clip(Omega_i_he, -wrn / Qn, 2.0e4))
        else:
            delta_w1_he = max(B1, 1.0)
            Omega_r_he = max(0.5 * dw_rn, 1.0)
            Omega_i_he = 500.0
            Om_init_he = complex(0.60 * dw_rn, 0.30 * dw_rn)
    else:
        B1 = (
            (dw_abs / 4.0)
            - (ws2 / (4.0 * dw_abs))
            + (xi * n_harmonic * Rn * wrn) / (16.0 * Qn * (dw_abs ** 2))
        )
        B2 = (xi * n_harmonic * wrn * Rn) / (4.0 * dw_abs * Qn)
        delta_w2_he = -B1 + np.sqrt(max(0.0, B1 * B1 + B2))
        if abs(delta_w2_he) < 1e-6:
            delta_w2_he = 1e-6
        delta_w1_he = -delta_w2_he
        Omega_r_he = dw_abs + delta_w2_he
        bn = -(xi * n_harmonic * Rn * (wrn ** 2)) / (4.0 * (Qn ** 2) * (delta_w2_he ** 2))
        kn = -(xi * n_harmonic * Rn * wrn) / (2.0 * Qn * (delta_w2_he ** 2))
        denom_he = 2.0 * Omega_r_he - kn - k1
        if abs(denom_he) < 1e-9:
            denom_he = -1e-9 if denom_he <= 0.0 else 1e-9
        Omega_i_he = (bn + b1 - (2.0 * Omega_r_he / tau_z_s)) / denom_he
        Om_init_he = complex(Omega_r_he, np.clip(Omega_i_he, -2.0 * wrn / Qn, 1.0e4))

    # 3. Guarded 20-step Newton refinement on exact coupled dispersion g(Omega) = 0
    gam_cav = wrn / (2.0 * Qn)
    if gam_cav < 0.25 * dw_abs:
        Om_D = Om_init_he
        max_step_D = (
            0.35 * dw_abs
            if disc < 0.0
            else max(0.35 * abs(delta_w1_he), 2.0 * gam_cav, 2.0 * np.pi * 80.0)
        )
    else:
        Om_D = complex(dw_abs, -gam_cav * 1.005)
        max_step_D = 0.25 * dw_abs

    for _ in range(20):
        g_val, dg_val = _eval_coupled_robinson_dispersion(
            Om_D, ws2, tau_z_s, xi, n_harmonic, R1_eff, Q1_eff, wr1, w10, Rn, Qn, wrn, wn0
        )
        if abs(dg_val) < 1e-24:
            break
        step = g_val / dg_val
        if abs(step) > max_step_D:
            step = step * (max_step_D / abs(step))
        Om_next = Om_D - step
        if abs(Om_next.real - dw_abs) > 0.65 * dw_abs or Om_next.real <= 0.0:
            break
        Om_D = Om_next
        if abs(step) < 1e-5:
            break

    # Also track C-mode (synchrotron-like dipole mode near sqrt(ws2_C))
    Om_C = complex(max(np.sqrt(ws0_sq), 2.0 * np.pi * 50.0), -1.0 / tau_z_s)
    for _ in range(20):
        g_c, dg_c = _eval_coupled_robinson_dispersion(
            Om_C, ws2_C, tau_z_s, xi, n_harmonic, R1_eff, Q1_eff, wr1, w10, Rn, Qn, wrn, wn0
        )
        if abs(dg_c) < 1e-24:
            break
        step_c = g_c / dg_c
        if abs(step_c) > 0.4 * abs(Om_C):
            step_c = step_c * (0.4 * abs(Om_C) / abs(step_c))
        Om_C = Om_C - step_c
        if abs(step_c) < 1e-5:
            break

    if (not np.isfinite(Om_D.real)) or abs(Om_D.real - dw_abs) > 0.60 * dw_abs:
        if gam_cav < 0.25 * dw_abs:
            Om_D = complex(max(Omega_r_he, 1.0), Omega_i_he)
        else:
            Om_D = complex(dw_abs, -gam_cav)

    Omega_r = float(abs(Om_D.real))
    Omega_i = float(Om_D.imag)
    delta_w1 = dw_abs - Omega_r

    return {
        "Omega_r": Omega_r,
        "Omega_i": Omega_i,
        "delta_w1": float(delta_w1),
        "Omega_r_he": float(Omega_r_he),
        "Omega_i_he": float(Omega_i_he),
        "Omega_C_r": float(Om_C.real),
        "Omega_C_i": float(Om_C.imag),
    }


def calculate_d_mode_robinson(
    I0: float,
    E_GeV: float,
    C_m: float,
    alpha_c: float,
    tau_z_s: float,
    U0_keV: float,
    Vrf_MV: float,
    R1_Q1: float,
    Q1: float,
    df_r1_Hz: float,
    n_harmonic: int,
    Rn_Qn: float,
    Qn: float,
    df_rn_Hz: float,
    h: Optional[int] = None,
    f_rf_Hz: Optional[float] = None,
    dfb_gain: float = 0.0,
    compute_threshold: bool = True,
) -> Dict[str, Any]:
    """Calculate D-mode Robinson instability frequency, growth rate, and critical detuning threshold.

    Implements T. He et al., Phys. Rev. Accel. Beams 27, 034402 (2024) [arXiv:2311.13205]
    with three physical enhancements:
    1. Explicit harmonic number ``h`` (fixing ``w0*C/(2*pi*c) == 1`` cancellation) and Eq. (5)
       equilibrium HHC synchrotron slope shift ``ws^2``.
    2. Closed-loop Direct RF Feedback (DFB) effective main-cavity impedance reduction ``(R1_eff, Q1_eff)``.
    3. Guarded Newton refinement on the exact coupled dipole dispersion relation ``g(Omega) = 0``,
       seamlessly supporting both high-Q superconducting PSHCs (e.g. HALF, Qn ~ 2e8) and
       low-Q warm copper HHCs (e.g. SOLEIL II, Qn ~ 3.1e4).
    """
    c = 2.99792458e8
    T0 = C_m / c
    f0 = 1.0 / T0
    w0 = 2.0 * np.pi * f0

    if h is not None:
        h_int = int(round(h))
    elif f_rf_Hz is not None and f_rf_Hz > 0.0:
        h_int = max(1, int(round(f_rf_Hz / f0)))
    elif abs(C_m - 353.98) < 2.0:
        h_int = 416
    else:
        h_int = max(1, int(round(499.8e6 / f0)))

    E0_eV = E_GeV * 1e9
    U0_eV = U0_keV * 1e3
    Vrf_V = Vrf_MV * 1e6

    R1 = float(R1_Q1) * float(Q1)
    Rn = float(Rn_Qn) * float(Qn)

    # Effective closed-loop MC impedance under Direct RF Feedback
    R1_eff, Q1_eff = R1, float(Q1)
    if dfb_gain > 0.0 and R1 > 0.0 and Q1 > 0.0:
        fr1 = h_int * f0 + df_r1_Hz
        psi1 = np.arctan(2.0 * Q1 * df_r1_Hz / max(fr1, 1.0))
        loop_factor = 1.0 + dfb_gain / max(np.cos(psi1), 0.1)
        R1_eff = R1 / loop_factor
        Q1_eff = float(Q1) / loop_factor

    xi = (h_int * max(I0, 1e-12) * alpha_c * w0) / (T0 * E0_eV)
    sin_phis = min(U0_eV / max(Vrf_V, 1.0), 0.9999)
    cos_phis_abs = np.sqrt(1.0 - sin_phis * sin_phis)
    ws0_sq = (xi * Vrf_V * cos_phis_abs) / max(I0, 1e-12)

    sol = _solve_dmode_at_detune(
        df_rn_Hz=df_rn_Hz,
        I0=I0,
        ws0_sq=ws0_sq,
        xi=xi,
        tau_z_s=tau_z_s,
        w0=w0,
        h=h_int,
        R1_eff=R1_eff,
        Q1_eff=Q1_eff,
        df_r1_Hz=df_r1_Hz,
        n_harmonic=n_harmonic,
        Rn=Rn,
        Qn=float(Qn),
    )

    # Optimal flat-potential (xi = 1.0) detuning (Hz)
    k_fp_sq = 1.0 / (n_harmonic ** 2) - (1.0 / max(n_harmonic ** 2 - 1.0, 1.0)) * ((U0_eV / max(Vrf_V, 1.0)) ** 2)
    V_hc_opt = np.sqrt(max(k_fp_sq, 0.0)) * Vrf_V
    frn0 = n_harmonic * h_int * f0
    if 2.0 * I0 * Rn > V_hc_opt and V_hc_opt > 0.0 and Qn > 0.0:
        cos_psi_opt = V_hc_opt / (2.0 * I0 * Rn)
        tan_psi_opt = np.sqrt(max(0.0, 1.0 / (cos_psi_opt * cos_psi_opt) - 1.0))
        opt_detune_Hz = float(frn0 * tan_psi_opt / (2.0 * Qn))
    else:
        opt_detune_Hz = 0.0

    # Critical D-mode threshold detuning where Omega_i(df_thresh) = 0
    thresh_detune_Hz = 0.0
    if compute_threshold and I0 > 1e-6 and Rn > 0.0 and Qn > 0.0:
        df_grid = np.logspace(2.0, 6.0, 48)  # 100 Hz to 1 MHz
        omi_grid = np.array(
            [
                _solve_dmode_at_detune(
                    df_g, I0, ws0_sq, xi, tau_z_s, w0, h_int,
                    R1_eff, Q1_eff, df_r1_Hz, n_harmonic, Rn, float(Qn)
                )["Omega_i"]
                for df_g in df_grid
            ]
        )
        pos_idx = np.where(omi_grid > 0.0)[0]
        if len(pos_idx) > 0:
            last_pos = int(pos_idx[-1])
            if last_pos < len(df_grid) - 1:
                a_df = float(df_grid[last_pos])
                b_df = float(df_grid[last_pos + 1])
                for _ in range(20):
                    m_df = 0.5 * (a_df + b_df)
                    omi_m = _solve_dmode_at_detune(
                        m_df, I0, ws0_sq, xi, tau_z_s, w0, h_int,
                        R1_eff, Q1_eff, df_r1_Hz, n_harmonic, Rn, float(Qn)
                    )["Omega_i"]
                    if omi_m > 0.0:
                        a_df = m_df
                    else:
                        b_df = m_df
                thresh_detune_Hz = 0.5 * (a_df + b_df)
            else:
                thresh_detune_Hz = float(df_grid[-1])

    Omega_i = sol["Omega_i"]
    return {
        "f_Dmode_Hz": sol["Omega_r"] / (2.0 * np.pi),
        "df1_Hz": sol["delta_w1"] / (2.0 * np.pi),
        "growth_rate_s1": Omega_i,
        "is_unstable": bool(Omega_i > 0.0),
        "f_Dmode_closed_form_Hz": sol["Omega_r_he"] / (2.0 * np.pi),
        "growth_rate_closed_form_s1": sol["Omega_i_he"],
        "f_Cmode_Hz": sol["Omega_C_r"] / (2.0 * np.pi),
        "cmode_growth_rate_s1": sol["Omega_C_i"],
        "thresh_detune_Hz": float(thresh_detune_Hz),
        "opt_detune_Hz": float(opt_detune_Hz),
        "is_stretching_safe": bool(opt_detune_Hz > thresh_detune_Hz),
    }


class InstabilityClassifier:
    """Comprehensive accelerator RF & beam instability classifier."""

    def __init__(self, params: Optional[MachineParameters] = None) -> None:
        self.params = params or MachineParameters()

    def solve_equilibrium_form_factor(
        self,
        I0: float,
        fre_shift_hc: float,
        use_form_factor: Optional[bool] = None,
    ) -> Dict[str, float]:
        """Self-consistent Fokker-Planck potential-well solver for stretched bunch form factor F_n(sigma_tau)."""
        p = self.params
        I0_A = I0 / 1000.0 if I0 > 5.0 else I0
        df_hz = fre_shift_hc * 1000.0 if abs(fre_shift_hc) < 1000.0 else fre_shift_hc
        use_ff = p.use_form_factor if use_form_factor is None else bool(use_form_factor)

        n = p.n_hc
        ratio = p.U0 / max(p.V_mc, 1.0)
        k_fp_sq = 1.0 / (n**2) - (1.0 / max(n**2 - 1.0, 1.0)) * (ratio**2)
        k_fp = float(np.sqrt(max(k_fp_sq, 1e-8)))

        fhc = p.fhc
        wre_hc = 2.0 * np.pi * fhc
        wr_hc = wre_hc + 2.0 * np.pi * df_hz
        psi_hc = float(np.arctan(p.Q_hc * (wr_hc / wre_hc - wre_hc / wr_hc)))

        Qs0 = max(p.synchrotron_tune_0, 1e-6)
        sig_t0 = p.alpha_c * p.sigma_delta0 / (2.0 * np.pi * Qs0) * p.T0

        if (not use_ff) or I0_A <= 1e-9 or p.R_hc <= 0 or p.Q_hc <= 0:
            xi_pt = (2.0 * I0_A * p.R_hc * abs(np.cos(psi_hc))) / max(k_fp * p.V_mc, 1.0)
            return {
                "F_hc": 1.0,
                "F_mc": 1.0,
                "sigma_t_eq": float(sig_t0),
                "sigma_t0": float(sig_t0),
                "xi_eff": float(xi_pt),
            }

        w_rf = 2.0 * np.pi * p.frf
        tau_max = min(0.45 / p.frf, max(12.0 * sig_t0, 320e-12))
        tau = np.linspace(-tau_max, tau_max, 301)
        dtau = tau[1] - tau[0]

        Fn_c = 0.90 + 0.0j
        F1_c = 0.994 + 0.0j
        sig_t = sig_t0
        converged = False
        phis_cent = np.pi - np.arcsin(min(0.999, p.U0 / max(p.V_mc, 1.0)))

        for it in range(30):
            arg_s = np.clip(
                (p.U0 + 2.0 * (abs(Fn_c) ** 2) * I0_A * p.R_hc * (np.cos(psi_hc) ** 2)) / max(abs(F1_c) * p.V_mc, 1.0),
                -0.999,
                0.999,
            )
            phis = np.pi - np.arcsin(arg_s) + float(np.angle(F1_c))
            E_gain = (
                p.V_mc * np.sin(phis + w_rf * tau)
                - 2.0 * I0_A * p.R_hc * np.cos(psi_hc) * np.real(np.exp(1j * psi_hc) * Fn_c * np.exp(1j * n * w_rf * tau))
                - p.U0
            )
            U_pot = -np.concatenate(([0.0], np.cumsum(0.5 * (E_gain[:-1] + E_gain[1:])) * dtau)) / (p.T0 * p.E0)
            U_pot -= np.min(U_pot)
            rho = np.exp(-np.clip(U_pot / max(p.alpha_c * (p.sigma_delta0**2), 1e-18), 0.0, 100.0))
            rho /= max(float(np.sum(rho) * dtau), 1e-30)

            t_mean = float(np.sum(tau * rho) * dtau)
            tau_c = tau - t_mean
            sig_t = float(np.sqrt(np.sum((tau_c**2) * rho) * dtau))

            F1_new = np.sum(rho * np.exp(-1j * w_rf * tau_c)) * dtau
            Fn_new = np.sum(rho * np.exp(-1j * n * w_rf * tau_c)) * dtau
            err = abs(Fn_new - Fn_c)
            Fn_c = 0.35 * Fn_new + 0.65 * Fn_c
            F1_c = 0.35 * F1_new + 0.65 * F1_c
            phis_cent = phis + w_rf * t_mean
            if err < 1e-6 and it >= 4:
                converged = True
                break

        xi_raw = (2.0 * abs(Fn_c) * I0_A * p.R_hc * abs(np.cos(psi_hc))) / max(k_fp * p.V_mc, 1.0)
        if (not converged) or phis_cent > (np.pi - 0.03) or phis_cent < (0.5 * np.pi) or (xi_raw > 1.06 and sig_t < 1.8 * sig_t0):
            # Over-stretched double-well PTBL regime: clamp form factor smoothly at flat-potential envelope (xi_target = 1.04)
            df_ref_khz = p.get_nominal_hc_detuning(I0_A * 1e3, target_xi=1.04, use_form_factor=True)
            psi_ref = np.arctan(2.0 * p.Q_hc * (df_ref_khz * 1e3) / fhc)
            Fn_mag = float(np.clip((1.04 * k_fp * p.V_mc) / max(2.0 * I0_A * p.R_hc * np.cos(psi_ref), 1e-9), 0.70, 1.0))
            F1_mag = 0.992
            sig_t = max(sig_t, 4.2 * sig_t0)
        else:
            Fn_mag = float(abs(Fn_c))
            F1_mag = float(abs(F1_c))

        xi_eff = (2.0 * Fn_mag * I0_A * p.R_hc * abs(np.cos(psi_hc))) / max(k_fp * p.V_mc, 1.0)
        return {
            "F_hc": Fn_mag,
            "F_mc": F1_mag,
            "sigma_t_eq": float(sig_t),
            "sigma_t0": float(sig_t0),
            "xi_eff": float(xi_eff),
        }

    def compute_flat_potential_parameter(
        self,
        I0: float,
        fre_shift_hc: float,
        use_form_factor: Optional[bool] = None,
    ) -> float:
        """Compute flat-potential distortion parameter xi (with optional iterative form factor F_n(sigma_tau)).

        xi = 2 * F_n * I0 * R_hc * cos(psi_hc) / (k_fp * V_mc)
        Automatically handles I0 in mA or A, and fre_shift_hc in kHz or Hz.
        """
        eq = self.solve_equilibrium_form_factor(I0, fre_shift_hc, use_form_factor=use_form_factor)
        return float(eq["xi_eff"])

    def analyze_llrf_nyquist(self, K_p: float, delay_sec: float) -> Tuple[bool, float, float]:
        """Perform physical LLRF frequency-domain open-loop Nyquist stability & margin evaluation.

        Evaluates the baseband open-loop transfer function matching MC_Cavity_Control.m & DFB_Set.m:
            L(j*w_m) = K_p * sec(psi_mc) / (1 + j*(w_m / w_half - tan(psi_mc))) * exp(-j*(w_m*delay_sec + phi_dfb))
        across both positive and negative modulation frequencies w_m, checking the Cauchy
        phase-crossover criterion (|L(j*w_pi)| < 1) and Euclidean vector margin min |1 + L(j*w_m)|.

        Returns:
            (is_nyquist_stable, normalized_margin, distance_to_minus_one)
        """
        if K_p < 0.0:
            return False, 0.0, 0.0
        if K_p == 0.0:
            return True, 1.0, 1.0

        p = self.params
        f_rf = max(p.frf, 1.0)
        q_mc = max(p.Q_mc, 1.0)
        w_half = 2.0 * np.pi * f_rf / (2.0 * q_mc)
        tan_psi = 2.0 * q_mc * p.fre_shift_mc / f_rf
        sec_psi = float(np.sqrt(1.0 + tan_psi * tan_psi))
        phi_dfb = float(p.DFB_phase_shift)
        tau_d = max(0.0, float(delay_sec))

        w_pos = np.geomspace(1e2, 5e7, 1200)
        L_pos = (K_p * sec_psi / (1.0 + 1j * (w_pos / w_half - tan_psi))) * np.exp(-1j * (w_pos * tau_d + phi_dfb))
        L_neg = (K_p * sec_psi / (1.0 + 1j * (-w_pos / w_half - tan_psi))) * np.exp(1j * (w_pos * tau_d + phi_dfb))

        ph_pos = np.unwrap(np.angle(L_pos))
        ph_neg = np.unwrap(np.angle(L_neg))

        gain_at_pi = 0.0
        idx_p = np.where(ph_pos <= -np.pi)[0]
        if len(idx_p) > 0:
            k = int(idx_p[0])
            if k > 0:
                frac = (-np.pi - ph_pos[k - 1]) / (ph_pos[k] - ph_pos[k - 1])
                gain_p = float(abs(L_pos[k - 1]) + frac * (abs(L_pos[k]) - abs(L_pos[k - 1])))
            else:
                gain_p = float(abs(L_pos[0]))
            gain_at_pi = max(gain_at_pi, gain_p)

        idx_n = np.where(ph_neg >= np.pi)[0]
        if len(idx_n) > 0:
            k = int(idx_n[0])
            if k > 0:
                frac = (np.pi - ph_neg[k - 1]) / (ph_neg[k] - ph_neg[k - 1])
                gain_n = float(abs(L_neg[k - 1]) + frac * (abs(L_neg[k]) - abs(L_neg[k - 1])))
            else:
                gain_n = float(abs(L_neg[0]))
            gain_at_pi = max(gain_at_pi, gain_n)

        dist_minus1 = float(min(np.min(np.abs(1.0 + L_pos)), np.min(np.abs(1.0 + L_neg))))
        is_stable = (gain_at_pi < 1.0) and (dist_minus1 > 0.05)
        normalized_margin = float(min(1.0, dist_minus1)) if is_stable else 0.0
        return bool(is_stable), normalized_margin, dist_minus1

    def evaluate_dmode_robinson(
        self,
        I0: float,
        fre_shift_hc: float,
        dfb_gain: Optional[float] = None,
        compute_threshold: bool = True,
    ) -> Dict[str, Any]:
        """Evaluate analytical D-mode Robinson frequency, growth rate, and threshold on ``self.params``."""
        p = self.params
        I0_A = I0 / 1000.0 if I0 > 5.0 else I0
        df_hz = fre_shift_hc * 1000.0 if abs(fre_shift_hc) < 1000.0 else fre_shift_hc
        g_dfb = dfb_gain if dfb_gain is not None else (p.DFB_gain if p.DFB_enable else 0.0)
        R1_Q1 = (p.R_mc / p.Q_mc) if p.Q_mc > 0.0 else p.RoQ_mc
        Rn_Qn = (p.R_hc / p.Q_hc) if p.Q_hc > 0.0 else p.RoQ_hc

        return calculate_d_mode_robinson(
            I0=I0_A,
            E_GeV=p.E0 * 1e-9,
            C_m=p.C,
            alpha_c=p.alpha_c,
            tau_z_s=p.tau_z,
            U0_keV=p.U0 * 1e-3,
            Vrf_MV=p.V_mc * 1e-6,
            R1_Q1=R1_Q1,
            Q1=p.Q_mc,
            df_r1_Hz=p.fre_shift_mc,
            n_harmonic=p.n_hc,
            Rn_Qn=Rn_Qn,
            Qn=p.Q_hc,
            df_rn_Hz=df_hz,
            h=p.h,
            dfb_gain=g_dfb,
            compute_threshold=compute_threshold,
        )

    def compute_dmode_threshold_curve(
        self,
        I0_mA_array: np.ndarray,
        dfb_gain: Optional[float] = None,
    ) -> Dict[str, np.ndarray]:
        """Scan critical D-mode threshold detuning (Omega_i = 0) vs optimal lengthening detuning across I0."""
        i0_mA = np.asarray(I0_mA_array, dtype=float)
        thresh_khz = np.zeros_like(i0_mA)
        opt_khz = np.zeros_like(i0_mA)
        is_safe = np.ones_like(i0_mA, dtype=bool)
        for idx, cur_mA in enumerate(i0_mA):
            dm = self.evaluate_dmode_robinson(
                cur_mA * 1e-3,
                max(self.params.fre_shift_hc, 10e3),
                dfb_gain=dfb_gain,
                compute_threshold=True,
            )
            thresh_khz[idx] = dm["thresh_detune_Hz"] * 1e-3
            opt_khz[idx] = dm["opt_detune_Hz"] * 1e-3
            is_safe[idx] = bool(dm["is_stretching_safe"])
        return {
            "I0_mA": i0_mA,
            "thresh_khz": thresh_khz,
            "opt_khz": opt_khz,
            "is_stretching_safe": is_safe,
        }


    def evaluate_from_signals(
        self,
        q_mean: np.ndarray,
        p_mean: np.ndarray,
        q_std: np.ndarray,
        p_std: np.ndarray,
        I0: float,
        fre_shift_hc: float,
        sigma_t0: Optional[float] = None,
        recor_step: int = 10,
        num_macroparticles: int = 10000,
    ) -> InstabilityResult:
        """Classify beam state from multi-bunch tracking time-series arrays.

        Args:
            q_mean: [n_turns, h] bunch centroid time deviation in sigma_t0 units
            p_mean: [n_turns, h] bunch centroid energy deviation
            q_std: [n_turns, h] bunch RMS length in sigma_t0 units
            p_std: [n_turns, h] bunch RMS energy spread
            I0: Beam current in [A]
            fre_shift_hc: Harmonic cavity detuning [Hz]
            sigma_t0: Natural bunch length [s]
            recor_step: Decimation turn step
            num_macroparticles: Macroparticle count per bunch
        """
        p = self.params
        sig_t0 = sigma_t0 or p.sigma_t0
        T0 = p.T0
        tau_z = p.tau_z
        damping_rate = p.damping_rate

        n_turns, h = q_mean.shape

        # 1. Check for extreme beam loss / NaN / Inf
        if np.any(np.isnan(q_mean)) or np.any(np.isinf(q_mean)) or np.max(np.abs(q_mean)) > 200.0:
            return InstabilityResult(
                status=StabilityStatus.LOST,
                is_stable=False,
                mode_id=0,
                growth_rate=1.0 / T0,
                bunch_length_factor=0.0,
                centroid_jitter_ps=999.0,
                energy_spread_factor=999.0,
                crest_factor=0.0,
                nyquist_margin=0.0,
                details={"reason": "Numerical divergence or catastrophic beam loss"},
            )

        # 2. Window selection
        idx_transient = int(n_turns * 0.5)
        idx_steady = int(n_turns * 0.8)
        win_analysis = slice(idx_transient, n_turns)
        win_steady = slice(idx_steady, n_turns)

        # 3. Bunch centroid jitter (ps)
        bunch_jitters = np.std(q_mean[win_steady, :], axis=0) * sig_t0 * 1e12
        centroid_jitter_ps = float(np.mean(bunch_jitters))

        # 4. Energy envelope and Hilbert growth rate
        P_detrend = p_mean[win_analysis, :] - np.mean(p_mean[win_analysis, :], axis=0)
        rms_P = np.sqrt(np.mean(P_detrend**2, axis=1))  # [turns_analysis]

        # Use Hilbert transform to extract smooth amplitude envelope
        analytic_signal = hilbert(rms_P - np.mean(rms_P))
        envelope = np.abs(analytic_signal) + np.mean(rms_P)

        turns_axis = np.arange(idx_transient, n_turns) * recor_step
        time_axis = turns_axis * T0

        noise_floor_P = 1.0 / np.sqrt(max(num_macroparticles, 100))
        noise_floor_ps = (sig_t0 * 1e12) / np.sqrt(max(num_macroparticles, 100))

        if np.max(envelope) > 2.0 * noise_floor_P:
            # Exponential growth rate fit on envelope
            log_env = np.log(np.maximum(envelope, 1e-15))
            poly = np.polyfit(time_axis, log_env, 1)
            growth_rate = float(poly[0])
        else:
            growth_rate = -0.5 * damping_rate

        # 5. AC-detrended Spatial Discrete Fourier Transform (DFT) for multi-bunch mode identification
        # Subtract static per-bunch equilibrium phase offset first so mu=0 reflects pure Mode-0 oscillation
        qm_steady_ac = q_mean[win_steady, :] - np.mean(q_mean[win_steady, :], axis=0, keepdims=True)
        if h > 1:
            spatial_fft = np.fft.fft(qm_steady_ac, axis=1) / float(h)
            mode_power = np.mean(np.abs(spatial_fft)**2, axis=0)  # length h
            mode0_pwr = float(mode_power[0])
            ac_modes = mode_power[1:]
            total_mb_pwr = float(np.sum(ac_modes))
            peak_mode_pwr = float(np.max(ac_modes))
            crest_factor = float(peak_mode_pwr / (np.mean(ac_modes) + 1e-30))
            peak_mb_rms_ps = float(np.sqrt(min(total_mb_pwr, 2.0 * peak_mode_pwr)) * sig_t0 * 1e12)
            is_multibunch_dom = (total_mb_pwr > mode0_pwr) and (peak_mb_rms_ps > max(1.5, 3.5 * noise_floor_ps))
            dom_mode_mu = int(np.argmax(ac_modes) + 1) if is_multibunch_dom else 0
        else:
            mode0_pwr = float(np.mean(qm_steady_ac**2))
            total_mb_pwr = 0.0
            crest_factor = 1.0
            peak_mb_rms_ps = 0.0
            is_multibunch_dom = False
            dom_mode_mu = 0

        has_lcbi_mode = is_multibunch_dom and (
            (crest_factor > 4.5 and (growth_rate > 0.3 * damping_rate or centroid_jitter_ps > 3.0))
            or (growth_rate > 0.5 * damping_rate and crest_factor > 3.5 and peak_mb_rms_ps > 2.0)
        )

        # 6. Bunch length stretching & oscillation
        bl_steady = q_std[win_steady, :] * sig_t0 * 1e12
        mean_bl_ps = float(np.mean(bl_steady))
        bunch_length_factor = mean_bl_ps / (sig_t0 * 1e12)

        mean_bl_t = np.mean(q_std[win_analysis, :], axis=1) * sig_t0 * 1e12
        bl_osc_rel = float(np.std(mean_bl_t) / (np.mean(mean_bl_t) + 1e-30))

        # 7. Energy spread factor
        mean_p_std = float(np.mean(p_std[win_steady, :]))
        energy_spread_factor = mean_p_std

        # 8. Physical potential parameter xi & Mode-0 coherent oscillation spectrum vs D-mode theory
        xi_calc = self.compute_flat_potential_parameter(I0, fre_shift_hc)

        sig_m0 = np.mean(q_mean[win_analysis, :], axis=1)
        sig_m0 = sig_m0 - np.mean(sig_m0)
        n_m0 = len(sig_m0)
        if n_m0 >= 8:
            win_hann = 0.5 * (1.0 - np.cos(2.0 * np.pi * np.arange(n_m0) / max(n_m0 - 1, 1)))
            n_pad = max(1024, 1 << int(np.ceil(np.log2(4 * n_m0))))
            sp_m0 = np.abs(np.fft.rfft(sig_m0 * win_hann, n=n_pad))
            f_ax0 = np.fft.rfftfreq(n_pad, d=recor_step * T0)
            f_mode0_Hz = float(f_ax0[int(np.argmax(sp_m0[1:]) + 1)]) if len(sp_m0) > 2 else 0.0
        else:
            f_mode0_Hz = 0.0

        dm_sig = self.evaluate_dmode_robinson(I0, fre_shift_hc, compute_threshold=False)

        # 9. Dynamic physical classification
        # Current-adaptive jitter threshold for flat-potential island
        jitter_tol_island = 5.0 + 10.0 * max(0.0, (I0 - 0.35) / 0.15)
        is_flat_stable = (
            (0.85 <= xi_calc <= 1.07)
            and (growth_rate <= 0.35 * damping_rate)
            and (energy_spread_factor < 1.25)
            and (bl_osc_rel < 0.08)
            and (centroid_jitter_ps < jitter_tol_island)
            and (not has_lcbi_mode)
        )

        if is_flat_stable:
            status = StabilityStatus.STABLE
        elif xi_calc >= 1.04 and (energy_spread_factor > 1.8 or bl_osc_rel > 0.15 or centroid_jitter_ps > 30.0):
            if I0 >= 0.365:
                status = StabilityStatus.PTBL
            else:
                status = StabilityStatus.COUPLED_DQ
        elif has_lcbi_mode:
            status = StabilityStatus.LCBI
        elif xi_calc >= 0.86 and (
            bl_osc_rel > 0.05
            or centroid_jitter_ps > max(3.5, 4.0 * noise_floor_ps)
            or growth_rate > 0.35 * damping_rate
        ):
            status = StabilityStatus.COUPLED_DQ
        elif (
            growth_rate > 0.35 * damping_rate
            or centroid_jitter_ps > max(5.0, 5.0 * noise_floor_ps)
            or energy_spread_factor > 1.25
            or bl_osc_rel > 0.08
        ):
            status = StabilityStatus.LOST  # Code 4: Mode-0 Robinson / LLRF loop instability
        else:
            status = StabilityStatus.STABLE

        is_stable = (status == StabilityStatus.STABLE)

        return InstabilityResult(
            status=status,
            is_stable=is_stable,
            mode_id=dom_mode_mu if not is_stable else 0,
            growth_rate=growth_rate,
            bunch_length_factor=bunch_length_factor,
            centroid_jitter_ps=centroid_jitter_ps,
            energy_spread_factor=energy_spread_factor,
            crest_factor=crest_factor,
            nyquist_margin=1.0 if is_stable else 0.0,
            details={
                "xi": xi_calc,
                "bl_osc_rel": bl_osc_rel,
                "mean_bl_ps": mean_bl_ps,
                "f_mode0_Hz": f_mode0_Hz,
                "f_Dmode_theory_Hz": dm_sig["f_Dmode_Hz"],
                "f_Cmode_theory_Hz": dm_sig["f_Cmode_Hz"],
                "dmode_growth_rate_theory_s1": dm_sig["growth_rate_s1"],
            },
        )

    def evaluate_semianalytical(
        self,
        I0: Optional[float] = None,
        fre_shift_hc: Optional[float] = None,
        *,
        I0_mA: Optional[float] = None,
        df_khz: Optional[float] = None,
        K_p: float = 0.25,
        tau_d: Optional[float] = None,
        delay_sec: Optional[float] = None,
    ) -> InstabilityResult:
        """Fast semi-analytical instability evaluation using Robinson equations and LLRF Nyquist.

        Supports both standard SI units (I0 [A], fre_shift_hc [Hz], tau_d [s])
        and accelerator user units (I0_mA [mA], df_khz [kHz], delay_sec [s]).
        Runs in sub-millisecond time.
        """
        p = self.params

        # Normalize beam current
        if I0_mA is not None:
            curr_mA = I0_mA
            curr_A = I0_mA / 1000.0
        elif I0 is not None:
            if I0 > 5.0:
                curr_mA = I0
                curr_A = I0 / 1000.0
            else:
                curr_A = I0
                curr_mA = I0 * 1000.0
        else:
            curr_A = p.I0
            curr_mA = p.I0 * 1000.0

        # Normalize detuning
        if df_khz is not None:
            detune_hz = df_khz * 1000.0
        elif fre_shift_hc is not None:
            detune_hz = fre_shift_hc * 1000.0 if abs(fre_shift_hc) < 1000.0 else fre_shift_hc
        else:
            detune_hz = p.fre_shift_hc

        # Normalize delay (defaults to machine's configured DFB loop delay)
        delay = delay_sec if delay_sec is not None else (tau_d if tau_d is not None else p.dfb_delay)

        eq_ff = self.solve_equilibrium_form_factor(curr_A, detune_hz)
        xi = float(eq_ff["xi_eff"])

        # 1. Nyquist check (when feedback loop is active, K_p > 0) & Analytical D-mode Robinson check
        if K_p > 0.0:
            is_nyq_stable, nyq_margin, _ = self.analyze_llrf_nyquist(K_p, delay)
        else:
            is_nyq_stable, nyq_margin = True, 1.0
        dm = self.evaluate_dmode_robinson(curr_A, detune_hz, dfb_gain=max(K_p, 0.0), compute_threshold=False)

        # 2. Analytical double-well / PTBL / D-mode Robinson check
        if xi > 1.05:
            if curr_A >= 0.365:
                status = StabilityStatus.PTBL
            else:
                status = StabilityStatus.COUPLED_DQ
            growth_rate = max(250.0, dm["growth_rate_s1"])
            is_stable = False
        elif dm["is_unstable"]:
            status = StabilityStatus.COUPLED_DQ if xi >= 0.85 else StabilityStatus.LOST
            growth_rate = float(dm["growth_rate_s1"])
            is_stable = False
        elif not is_nyq_stable:
            status = StabilityStatus.COUPLED_DQ
            growth_rate = 150.0
            is_stable = False
        elif 0.85 <= xi <= 1.05:
            status = StabilityStatus.STABLE
            growth_rate = max(-0.5 * p.damping_rate, float(dm["growth_rate_s1"]))
            is_stable = True
        elif xi < 0.85:
            status = StabilityStatus.STABLE
            growth_rate = max(-p.damping_rate, float(dm["growth_rate_s1"]))
            is_stable = True
        else:
            status = StabilityStatus.COUPLED_DQ
            growth_rate = 50.0
            is_stable = False

        # Bunch lengthening factor from self-consistent Fokker-Planck potential well (when enabled)
        if p.use_form_factor and eq_ff["sigma_t0"] > 0.0 and xi <= 1.05:
            bl_factor = max(1.0, eq_ff["sigma_t_eq"] / eq_ff["sigma_t0"])
        elif xi <= 1.0:
            bl_factor = 1.0 + 3.2 * (xi**2)
        else:
            bl_factor = max(1.0, 4.2 - 2.5 * (xi - 1.0))

        return InstabilityResult(
            status=status,
            is_stable=is_stable,
            mode_id=1 if status == StabilityStatus.PTBL else 0,
            growth_rate=growth_rate,
            bunch_length_factor=float(bl_factor),
            centroid_jitter_ps=3.0 if is_stable else 18.0,
            energy_spread_factor=1.0 if is_stable else 1.6,
            crest_factor=1.5 if is_stable else 6.0,
            nyquist_margin=nyq_margin,
            details={
                "xi": xi,
                "f_Dmode_Hz": dm["f_Dmode_Hz"],
                "df1_Hz": dm["df1_Hz"],
                "dmode_growth_rate_s1": dm["growth_rate_s1"],
            },
        )
