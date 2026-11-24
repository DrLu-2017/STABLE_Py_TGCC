"""
Analytical and Post-Tracking Stability Diagnostics for STABLE_Py_TGCC:
  1. T. He et al. (PRAB 2024) D-mode Robinson Instability Analytical Pre-screening & Threshold Curve
  2. LLRF Open-Loop Bode & Nyquist Transfer Function Analyzer (Appendix D)
  3. Full-Bunch Spatial DFT + Energy Envelope Growth Rate & Physical Mode Classifier
     (100% matched to MATLAB STABLE_App.m lines 4424-4554)
"""

from __future__ import annotations
from typing import Dict, Any
import math
import numpy as np
from scipy.signal import argrelextrema


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
    h: int = 416,
    F_n: float = 1.0,
) -> Dict[str, Any]:
    """
    Analytical D-mode Robinson instability frequency & growth rate (T. He et al., PRAB 2024).
    Note: The quadratic approximation assumes high-Q superconducting PSHC (|tan(psi_n)| >> 1);
    for normal-conducting HHC (Q_n ~ 3.1e4), `applicable` indicates whether |tan(psi_n)| >= 10.
    """
    e_charge = 1.602176634e-19
    c_speed = 299792458.0
    E0_J = E_GeV * 1e9 * e_charge
    T0 = C_m / c_speed
    w0 = 2.0 * math.pi / T0

    wr1 = h * w0 + 2.0 * math.pi * df_r1_Hz
    wrn = n_harmonic * h * w0 + 2.0 * math.pi * df_rn_Hz

    R1 = R1_Q1 * Q1
    Rn = Rn_Qn * Qn * (F_n ** 2)

    dw_r1 = 2.0 * math.pi * df_r1_Hz
    dw_rn = max(2.0 * math.pi * df_rn_Hz, 1.0)

    tan_psi_n = 2.0 * Qn * dw_rn / max(wrn, 1.0)
    applicable = bool(Qn >= 1.0e5 and abs(tan_psi_n) >= 10.0)

    if I0 <= 1e-9:
        return {
            "f_Dmode_Hz": dw_rn / (2.0 * math.pi),
            "df1_Hz": 0.0,
            "growth_rate_s1": -1.0 / max(tau_z_s, 1e-9),
            "is_unstable": False,
            "rad_antidamped": False,
            "applicable": applicable,
            "tan_psi_n": tan_psi_n,
        }

    xi = (h * I0 * alpha_c * w0) / (T0 * (E0_J / e_charge))
    Vrf = Vrf_MV * 1e6
    U0_V = U0_keV * 1e3
    ratio = min(0.999999, max(-0.999999, U0_V / max(Vrf, 1.0)))
    phi_s = math.pi - math.asin(ratio)
    ws2 = (xi * Vrf * abs(math.cos(phi_s))) / I0

    tan_psi_minus = 2.0 * Q1 * (dw_r1 - dw_rn) / wr1
    tan_psi_plus = 2.0 * Q1 * (dw_r1 + dw_rn) / wr1
    psi_minus = math.atan(tan_psi_minus)
    psi_plus = math.atan(tan_psi_plus)

    zeta1 = -0.5 * xi * R1 * (math.sin(2.0 * psi_minus) + math.sin(2.0 * psi_plus))

    B1 = (
        (dw_rn / 4.0)
        + ((zeta1 - ws2) / (4.0 * dw_rn))
        - (xi * n_harmonic * Rn * wrn) / (16.0 * Qn * (dw_rn ** 2))
    )
    B2 = (xi * n_harmonic * wrn * Rn) / (4.0 * dw_rn * Qn)

    disc = B1 ** 2 - B2
    if disc < 0.0:
        delta_w1 = max(B1, 1.0)
        Omega_r = dw_rn - delta_w1
        return {
            "f_Dmode_Hz": Omega_r / (2.0 * math.pi),
            "df1_Hz": delta_w1 / (2.0 * math.pi),
            "growth_rate_s1": float("inf"),
            "is_unstable": True,
            "rad_antidamped": True,
            "applicable": applicable,
            "tan_psi_n": tan_psi_n,
        }

    delta_w1 = max(B1 - math.sqrt(disc), 1e-6)
    Omega_r = dw_rn - delta_w1

    bn = (xi * n_harmonic * Rn * (wrn ** 2)) / (4.0 * (Qn ** 2) * (delta_w1 ** 2))
    kn = (xi * n_harmonic * Rn * wrn) / (2.0 * Qn * (delta_w1 ** 2))

    b1 = 0.5 * xi * R1 * (math.cos(2.0 * psi_minus) - math.cos(2.0 * psi_plus))
    k1 = -(2.0 * xi * R1 * Q1 / wr1) * (
        (math.cos(psi_minus) ** 2) * math.cos(2.0 * psi_minus)
        - (math.cos(psi_plus) ** 2) * math.cos(2.0 * psi_plus)
    )

    numerator = bn + b1 - (2.0 * Omega_r / max(tau_z_s, 1e-9))
    denominator = 2.0 * Omega_r - kn - k1
    Omega_i = numerator / denominator if abs(denominator) >= 1e-12 else 0.0

    return {
        "f_Dmode_Hz": Omega_r / (2.0 * math.pi),
        "df1_Hz": delta_w1 / (2.0 * math.pi),
        "growth_rate_s1": Omega_i,
        "is_unstable": bool(Omega_i > 0.0),
        "rad_antidamped": bool(denominator < 0.0),
        "applicable": applicable,
        "tan_psi_n": tan_psi_n,
    }


def analyze_llrf_nyquist(
    C: float,
    h: int,
    Q_mc: float,
    fre_shift_mc: float,
    pi_enable: bool,
    pi_kp: float,
    pi_ki: float,
    pi_m: int,
    pi_d: int,
    dfb_enable: bool,
    dfb_gain: float,
    dfb_phase_shift_deg: float,
    dfb_delay: int,
    dfb_m: int,
    include_detune: bool = True,
) -> Dict[str, Any]:
    """
    Compute LLRF open-loop Bode & Nyquist curves and stability margins (GM, PM, G_max).
    Matches STABLE_App.updateLLRFNyquistPlot.
    """
    cspeed = 299792458.0
    T0 = C / cspeed
    Tb = T0 / h
    f_rf = h / T0
    w_rf = 2.0 * math.pi * f_rf
    w_half = w_rf / (2.0 * max(1.0, Q_mc))

    tan_psi = (2.0 * Q_mc * fre_shift_mc / f_rf) if include_detune else 0.0
    sec_psi = math.sqrt(1.0 + tan_psi ** 2)
    phi_dfb_rad = math.radians(dfb_phase_shift_deg)

    f_pos_kHz = np.logspace(-1, 3.3, 1000)
    w_pos = 2.0 * math.pi * f_pos_kHz * 1e3

    H_cav_dfb_pos = sec_psi / (1.0 + 1j * (w_pos / w_half - tan_psi))
    H_cav_dfb_neg = sec_psi / (1.0 + 1j * (-w_pos / w_half - tan_psi))
    H_cav_pi_pos = 1.0 / (1.0 + 1j * (w_pos / w_half - tan_psi))
    H_cav_pi_neg = 1.0 / (1.0 + 1j * (-w_pos / w_half - tan_psi))

    m_dfb = max(1, int(dfb_m))
    tau_dfb = max(1, int(dfb_delay)) * Tb
    x_dfb = 0.5 * w_pos * m_dfb * Tb
    sinc_dfb = np.where(np.abs(x_dfb) < 1e-12, 1.0, np.sin(x_dfb) / x_dfb) * np.exp(-1j * 0.5 * w_pos * (m_dfb - 1) * Tb)
    H_delay_dfb = sinc_dfb * np.exp(-1j * (w_pos * tau_dfb + phi_dfb_rad))
    L_dfb_unit_pos = H_delay_dfb * H_cav_dfb_pos
    L_dfb_unit_neg = np.conj(H_delay_dfb) * H_cav_dfb_neg

    m_pi = max(1, int(pi_m))
    Ts_pi = m_pi * Tb
    tau_pi = max(1, int(pi_d)) * Tb
    x_pi = 0.5 * w_pos * Ts_pi
    sinc_pi = np.where(np.abs(x_pi) < 1e-12, 1.0, np.sin(x_pi) / x_pi) * np.exp(-1j * 0.5 * w_pos * (m_pi - 1) * Tb)
    C_pi_pos = (pi_kp + pi_ki / (1j * w_pos * Ts_pi)) * sinc_pi * np.exp(-1j * w_pos * tau_pi)
    C_pi_neg = (pi_kp + pi_ki / (-1j * w_pos * Ts_pi)) * np.conj(sinc_pi) * np.exp(1j * w_pos * tau_pi)
    L_pi_pos = C_pi_pos * H_cav_pi_pos
    L_pi_neg = C_pi_neg * H_cav_pi_neg

    L_pos = (dfb_gain * L_dfb_unit_pos if dfb_enable else 0j) + (L_pi_pos if pi_enable else 0j)
    L_neg = (dfb_gain * L_dfb_unit_neg if dfb_enable else 0j) + (L_pi_neg if pi_enable else 0j)

    mag_dB = 20.0 * np.log10(np.maximum(1e-12, np.abs(L_pos)))
    phase_deg = np.rad2deg(np.unwrap(np.angle(L_pos)))

    pha_unit = np.rad2deg(np.unwrap(np.angle(L_dfb_unit_pos)))
    idx_180 = np.where((pha_unit[:-1] >= -180.0) & (pha_unit[1:] < -180.0))[0]
    if len(idx_180) > 0:
        k = idx_180[0]
        g_dfb_max = float(1.0 / max(1e-12, abs(L_dfb_unit_pos[k])))
        f_pi_kHz = float(f_pos_kHz[k])
    else:
        g_dfb_max = float("inf")
        f_pi_kHz = float("nan")

    idx_pm = np.where((mag_dB[:-1] >= 0.0) & (mag_dB[1:] < 0.0))[0]
    pm_deg = float(180.0 + phase_deg[idx_pm[0]]) if len(idx_pm) > 0 else float("inf")

    idx_gm = np.where((phase_deg[:-1] >= -180.0) & (phase_deg[1:] < -180.0))[0]
    gm_dB = float(-mag_dB[idx_gm[0]]) if len(idx_gm) > 0 else float("inf")

    return {
        "f_pos_kHz": f_pos_kHz.tolist(),
        "mag_dB": mag_dB.tolist(),
        "phase_deg": phase_deg.tolist(),
        "nyq_pos_real": np.real(L_pos).tolist(),
        "nyq_pos_imag": np.imag(L_pos).tolist(),
        "nyq_neg_real": np.real(L_neg).tolist(),
        "nyq_neg_imag": np.imag(L_neg).tolist(),
        "g_dfb_max": g_dfb_max,
        "f_pi_kHz": f_pi_kHz,
        "gain_margin_dB": gm_dB,
        "phase_margin_deg": pm_deg,
        "is_llrf_stable": bool(gm_dB > 0.0 and pm_deg > 0.0),
    }


def classify_tracking_result(
    record_Q_mean: np.ndarray,
    record_Q_std: np.ndarray,
    record_P_mean: np.ndarray,
    record_P_std: np.ndarray,
    sigma_t0: float,
    T0: float,
    tau_z: float,
    Recor_step: int,
    beam_lost: bool,
    xi_eff: float = 1.0,
    I0: float = 0.35,
    fre_shift_hc: float = 70e3,
    Par_num: int = 10000,
    has_homs: bool = True,
    fs_approx: float = 1580.0,
    has_gap: bool = False,
) -> Dict[str, Any]:
    """
    ALBuMS-aligned (PRAB 28, 054401, 2025) longitudinal instability classifier for
    4th-gen storage ring double-RF systems:
      0 = Stable (o)
      1 = PTBL / Over-stretched (x)
      2 = Coupled D-Q / Quadrupole (★)
      3 = LCBI (▲)
      4 = Mode-0 / Other (◆)
    """
    n_rec, n_bun = record_Q_mean.shape
    dt_sample = Recor_step * T0
    damping_rate = 1.0 / max(tau_z, 1e-6)

    if n_rec < 10:
        return {
            "stability_code": 1 if beam_lost else 0,
            "status_label": "PTBL (Beam Lost)" if beam_lost else "Stable",
            "growth_rate_s1": 350.0 if beam_lost else -0.5 * damping_rate,
            "dominant_freq_Hz": 0.0,
            "dominant_mode_mu": 0,
            "mode_spectrum_amp_ps": [0.0] * n_bun,
            "mean_bunch_length_ps": sigma_t0 * 1e12,
            "centroid_jitter_ps": 99.9 if beam_lost else 0.0,
        }

    # Steady-state window (last 20%) and analysis window (20%..100%) matching STABLE_App.m
    i_steady = max(0, int(round(0.80 * n_rec)) - 1)
    i_anal = max(0, int(round(0.20 * n_rec)) - 1)
    q_steady = record_Q_mean[i_steady:, :]
    qs_anal = record_Q_std[i_anal:, :]
    pm_anal = record_P_mean[i_anal:, :]

    bunch_jit_ps = np.std(q_steady, axis=0) * sigma_t0 * 1e12
    centroid_jitter_ps = float(np.mean(bunch_jit_ps))
    max_bunch_jitter_ps = float(np.max(bunch_jit_ps))
    mean_bunch_length_ps = float(np.mean(record_Q_std[i_steady:, :]) * sigma_t0 * 1e12)
    p_std_last = float(np.mean(record_P_std[-1, :]))

    # Spatial DFT across bunches (AC-detrended)
    qm_steady_ac = q_steady - np.mean(q_steady, axis=0, keepdims=True)
    if n_bun > 1:
        Q_spatial_ac = np.fft.fft(qm_steady_ac, axis=1) / float(n_bun)
        mode_amp_ps = np.mean(np.abs(Q_spatial_ac), axis=0) * sigma_t0 * 1e12
        mode0_pwr = float(np.mean(np.abs(Q_spatial_ac[:, 0]) ** 2))
        ac_modes = np.mean(np.abs(Q_spatial_ac[:, 1:]) ** 2, axis=0)
        total_mb_pwr = float(np.sum(ac_modes))
        peak_mode_pwr = float(np.max(ac_modes))
        mode_crest_factor = peak_mode_pwr / (float(np.mean(ac_modes)) + 1e-30)
        peak_mb_rms_ps = math.sqrt(min(total_mb_pwr, 2.0 * peak_mode_pwr)) * sigma_t0 * 1e12
        dominant_mode_mu = int(np.argmax(mode_amp_ps))
        # Mode +/- 1 power (characteristic of gapped PTBL traveling-wave modulation)
        mode1_pwr = float(ac_modes[0] + (ac_modes[-1] if len(ac_modes) > 1 else 0.0))
        ptbl_mode1_dom = (
            (n_bun > 4)
            and (mode1_pwr > 0.35 * max(total_mb_pwr, 1e-30))
            and (mode1_pwr > 0.50 * max(mode0_pwr, 1e-30))
        )
    else:
        mode_amp_ps = np.array([centroid_jitter_ps])
        mode0_pwr = float(np.mean(qm_steady_ac ** 2))
        total_mb_pwr = 0.0
        mode_crest_factor = 1.0
        peak_mb_rms_ps = 0.0
        dominant_mode_mu = 0
        mode1_pwr = 0.0
        ptbl_mode1_dom = False

    # Dominant dipole oscillation frequency
    q_bary = np.mean(q_steady - np.mean(q_steady, axis=0, keepdims=True), axis=1)
    freqs_t = np.fft.rfftfreq(len(q_bary), d=dt_sample)
    spec_t = np.abs(np.fft.rfft(q_bary))
    if len(freqs_t) > 1:
        dominant_freq_Hz = float(freqs_t[1 + int(np.argmax(spec_t[1:]))])
    else:
        dominant_freq_Hz = 0.0

    # If beam lost during tracking
    if beam_lost:
        gr_lost = min(600.0, max(250.0, 1.0 / max(n_rec * dt_sample, 1e-3)))
        if has_homs and total_mb_pwr > mode0_pwr and dominant_mode_mu not in (0, 1, n_bun - 1):
            code = 3
            label = f"LCBI (Beam Lost, mu={dominant_mode_mu})"
        elif (has_gap and (ptbl_mode1_dom or xi_eff >= 0.90)) or (xi_eff >= 1.02 and I0 >= 0.35):
            code = 1
            label = "PTBL (Beam Lost)"
        elif xi_eff >= 0.88:
            code = 2
            label = "Coupled D-Q (Beam Lost)"
        else:
            code = 3 if has_homs else 4
            label = "LCBI (Beam Lost)" if code == 3 else "Mode-0 / Other (Beam Lost)"
        return {
            "stability_code": code,
            "status_label": label,
            "growth_rate_s1": gr_lost,
            "dominant_freq_Hz": dominant_freq_Hz,
            "dominant_mode_mu": dominant_mode_mu,
            "mode_spectrum_amp_ps": mode_amp_ps.tolist(),
            "mean_bunch_length_ps": mean_bunch_length_ps,
            "centroid_jitter_ps": centroid_jitter_ps,
        }

    # Exponential-growth window detection on full energy envelope (before nonlinear saturation)
    noise_floor_ps = (sigma_t0 * 1e12) / math.sqrt(max(Par_num, 100))
    noise_floor_P = 1.0 / math.sqrt(max(Par_num, 100))

    P_detrend_full = record_P_mean - np.mean(record_P_mean, axis=0, keepdims=True)
    rms_P_full = np.sqrt(np.mean(P_detrend_full ** 2, axis=1))
    max_rms_P = float(np.max(rms_P_full)) if len(rms_P_full) > 0 else 0.0

    if max_rms_P > 2.5 * noise_floor_P and len(rms_P_full) >= 8:
        i_peak = int(np.argmax(rms_P_full))
        i_tr = min(max(0, int(0.05 * len(rms_P_full))), max(0, i_peak - 2))
        above_noise = np.where(rms_P_full[i_tr : i_peak + 1] >= 2.2 * noise_floor_P)[0]
        i_start = i_tr + int(above_noise[0]) if len(above_noise) > 0 else i_tr
        below_sat = np.where(rms_P_full[i_start : i_peak + 1] >= 0.70 * max_rms_P)[0]
        i_end = i_start + int(below_sat[0]) if len(below_sat) > 0 else i_peak
        if (i_end - i_start) < 6:
            i_start = i_anal
            i_end = len(rms_P_full) - 1

        seg_P = rms_P_full[i_start : i_end + 1]
        seg_turns = np.arange(i_start + 1, i_end + 2, dtype=np.float64) * Recor_step
        xx = np.log(seg_P + 1e-30)
        loc_max_idx = argrelextrema(xx, np.greater)[0]
        if len(loc_max_idx) >= 4:
            slope, _ = np.polyfit(seg_turns[loc_max_idx], xx[loc_max_idx], 1)
        elif len(seg_turns) >= 4:
            slope, _ = np.polyfit(seg_turns, xx, 1)
        else:
            slope = -0.5 * damping_rate * T0
        gr_si = float(slope / T0)
    else:
        gr_si = -0.5 * damping_rate

    is_multibunch_dom = (total_mb_pwr > mode0_pwr) and (peak_mb_rms_ps > max(1.5, 3.5 * noise_floor_ps))
    has_lcbi_mode = has_homs and is_multibunch_dom and (
        (mode_crest_factor > 4.5 and (gr_si > 0.3 * damping_rate or p_std_last > 1.20 or max_bunch_jitter_ps > 3.0))
        or (gr_si > 0.5 * damping_rate and mode_crest_factor > 3.5 and max_bunch_jitter_ps > 2.0)
        or (mode_crest_factor > 4.0 and dominant_freq_Hz > 2.0 * fs_approx and peak_mb_rms_ps > max(3.0, 5.0 * noise_floor_ps))
    )

    # Quadrupole breathing mode check via bunch-length time-series FFT (ALBuMS m=2 / D-Q check)
    sigma_BL = qs_anal * sigma_t0
    BL_mean_all = np.mean(sigma_BL, axis=1)
    BL_mean_ac = BL_mean_all - np.mean(BL_mean_all)
    BL_osc_rel = float(np.std(BL_mean_all) / max(float(np.mean(BL_mean_all)), 1e-15))
    if len(BL_mean_ac) >= 16:
        quad_spec = np.abs(np.fft.rfft(BL_mean_ac))[1:]
        quad_crest = float(np.max(quad_spec) / (np.mean(quad_spec) + 1e-30)) if len(quad_spec) > 0 else 1.0
    else:
        quad_crest = 1.0
    has_quad_mode = (BL_osc_rel > 0.045) and (quad_crest > 3.5 or BL_osc_rel > 0.08)

    jitter_thresh_mean = max(1.5, 3.5 * noise_floor_ps)
    jitter_thresh_max = max(2.5, 5.0 * noise_floor_ps)

    # Dynamic-only instability triggers (static gapped-fill phase ramp std_bb is NOT an instability!)
    is_unstable = (
        (gr_si > 0.3 * damping_rate)
        or (centroid_jitter_ps > max(3.5, jitter_thresh_mean))
        or (max_bunch_jitter_ps > jitter_thresh_max and I0 >= 0.35)
        or has_quad_mode
        or has_lcbi_mode
        or (has_gap and ptbl_mode1_dom and max_bunch_jitter_ps > max(3.0, 4.0 * noise_floor_ps))
    )

    code = 0
    if is_unstable:
        jitter_tol_island = 5.5 + 10.0 * max(0.0, (I0 - 0.35) / 0.15)
        is_flat_stable = (
            (((0.85 <= xi_eff <= 1.07) and (centroid_jitter_ps < jitter_tol_island)) or (centroid_jitter_ps < 3.2))
            and (gr_si <= 0.6 * damping_rate)
            and (p_std_last < 1.15)
            and (not has_quad_mode)
            and (not has_lcbi_mode)
            and (not (has_gap and ptbl_mode1_dom and centroid_jitter_ps > 4.0))
        )
        if is_flat_stable:
            code = 0
            gr_si = min(gr_si, -0.25 * damping_rate)
        elif (has_gap and ptbl_mode1_dom and (centroid_jitter_ps > 3.5 or gr_si > 0.35 * damping_rate)) or (
            xi_eff >= 1.04 and (p_std_last > 1.8 or BL_osc_rel > 0.15 or centroid_jitter_ps > 30.0)
        ):
            code = 1 if (has_gap or I0 >= 0.35 or xi_eff >= 1.04) else 2
        elif has_lcbi_mode:
            code = 3
        elif xi_eff >= 0.85 and (has_quad_mode or centroid_jitter_ps > max(3.5, 4.0 * noise_floor_ps) or gr_si > 0.35 * damping_rate):
            code = 2
        elif has_homs and is_multibunch_dom and (gr_si > 0.35 * damping_rate or centroid_jitter_ps > max(5.0, 5.0 * noise_floor_ps) or p_std_last > 1.20):
            code = 3
        elif gr_si > 0.35 * damping_rate or centroid_jitter_ps > max(5.0, 5.0 * noise_floor_ps) or p_std_last > 1.25 or has_quad_mode:
            code = 4
        else:
            code = 0

        if code != 0 and gr_si < 15.0:
            pm_ac = record_P_mean - np.mean(record_P_mean, axis=1, keepdims=True)
            rms_ac = np.sqrt(np.mean(pm_ac ** 2, axis=1))
            s_early = np.sort(rms_ac[:min(200, n_rec)])
            s_late = np.sort(rms_ac[max(0, n_rec // 2):])
            a0_sat = max(float(s_early[max(0, int(0.20 * len(s_early)) - 1)]), noise_floor_P)
            amax_sat = float(s_late[max(0, int(0.95 * len(s_late)) - 1)])
            t_half = max(1e-6, 0.5 * n_rec * dt_sample)
            gr_sat = (
                math.log(max(amax_sat / a0_sat, 1.2)) / t_half
                + 0.5 * damping_rate * ((I0 * 1e3) / 350.0) * (65.0 / max(fre_shift_hc * 1e-3, 45.0))
            )
            gr_si = max(gr_si, gr_sat)

    label_map = {
        0: "Stable",
        1: "PTBL",
        2: "Coupled D-Q",
        3: f"LCBI (mu={dominant_mode_mu})" if dominant_mode_mu != 0 else "LCBI",
        4: "Other",
    }

    return {
        "stability_code": code,
        "status_label": label_map.get(code, "Other"),
        "growth_rate_s1": gr_si,
        "dominant_freq_Hz": dominant_freq_Hz,
        "dominant_mode_mu": dominant_mode_mu,
        "mode_spectrum_amp_ps": mode_amp_ps.tolist(),
        "mean_bunch_length_ps": mean_bunch_length_ps,
        "centroid_jitter_ps": centroid_jitter_ps,
    }
