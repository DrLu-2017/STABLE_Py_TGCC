"""
Post-processing module for extracting coupled-bunch instability mode and intrinsic
growth rate from multi-bunch tracking data.

Reference:
    "Beam instability induced by the rf cavity accelerating modes in STCF", Appendix C.
"""

from typing import Tuple, Optional
import numpy as np
import scipy.signal


def extract_tracking_instability(
    q: np.ndarray,
    tau_rad: float,
    downsample_interval: float = 10.0,
    fit_range: Optional[Tuple[int, int]] = None,
) -> Tuple[int, float, float, np.ndarray]:
    """
    Extract dominant coupled-bunch instability mode number and intrinsic growth rate
    from longitudinal bunch centroid tracking data (Appendix C algorithm).

    Parameters
    ----------
    q : np.ndarray, shape (N_samples, M_bunches)
        Longitudinal centroid positions q[n, m] for M bunches over N recorded samples.
    tau_rad : float
        Longitudinal radiation damping time (in the same time units as `downsample_interval`,
        e.g., turns if `downsample_interval` is in turns, or seconds if `downsample_interval`
        is `Recor_step * T0`).
    downsample_interval : float, optional
        Sampling interval between consecutive rows of `q` (default: 10 turns).
    fit_range : tuple of (int, int), optional
        Optional `(i_start, i_end)` slice indices within `A_envelope` specifying the
        exponential (log-linear) growth stage. If None, automatically selects the valid
        unsaturated log-linear segment of `A_envelope`.

    Returns
    -------
    mu_star : int
        Dominant coupled-bunch mode index in `[0, M - 1]` (`+mu_cand` for positive mode,
        or `M - mu_cand` corresponding to `-mu_cand` for negative mode).
    net_growth_rate : float
        Intrinsic (radiation-damping-compensated) growth rate:
        `net_growth_rate = growth_fit + 1.0 / tau_rad`.
    growth_fit : float
        Fitted exponential growth rate `1 / tau_fit` from `A_envelope`.
    A_envelope : np.ndarray, shape (N_window,)
        Time-domain amplitude envelope of the dominant mode over the last 1/3 window.
    """
    q = np.asarray(q, dtype=np.float64)
    if q.ndim != 2:
        raise ValueError(f"Expected 2D array (N_samples, M_bunches), got shape {q.shape}.")
    if tau_rad <= 0:
        raise ValueError(f"tau_rad must be positive, got {tau_rad}.")

    N, M = q.shape
    if N < 6 or M < 2:
        raise ValueError(f"Insufficient dimensions (N={N}, M={M}) for instability extraction.")

    # -------------------------------------------------------------------------
    # Step 1: Data Windowing (last 1/3 of time samples)
    # -------------------------------------------------------------------------
    start_idx = int(2 * N // 3)
    q_window = q[start_idx:, :]
    N_win = q_window.shape[0]

    # -------------------------------------------------------------------------
    # Step 2: Subtract Mode-0 (DC / barycenter) Component at Final Sample
    # -------------------------------------------------------------------------
    q_last = q_window[-1, :]
    q_m = q_last - np.mean(q_last)

    # -------------------------------------------------------------------------
    # Step 3: Spatial DFT to Identify Candidate Mode |mu_cand|
    # -------------------------------------------------------------------------
    dft_spectrum = np.abs(np.fft.fft(q_m)) / M
    # Explicitly zero out DC (index 0) to guarantee DC exclusion in argmax
    half_spectrum = dft_spectrum[: M // 2].copy()
    half_spectrum[0] = 0.0
    mu_cand = int(np.argmax(half_spectrum))

    # -------------------------------------------------------------------------
    # Step 4: Hilbert Transform & Complex Analytic Signal Construction
    # -------------------------------------------------------------------------
    # Note: scipy.signal.hilbert(x, axis=0) returns x + 1j * H[x], where H[x] is
    # the Hilbert transform of x. Therefore:
    #   q_analytic = q_window - 1j * H[q_window] = conj(scipy.signal.hilbert(q_window, axis=0))
    hilbert_transform_q = np.imag(scipy.signal.hilbert(q_window, axis=0))
    q_analytic = q_window - 1j * hilbert_transform_q

    # -------------------------------------------------------------------------
    # Step 5: Spatial Mode Projection & Positive/Negative Mode Disambiguation
    # -------------------------------------------------------------------------
    m_idx = np.arange(M, dtype=np.float64)[np.newaxis, :]  # shape (1, M)
    spatial_pos = np.exp(-1j * 2.0 * np.pi * (+mu_cand) * m_idx / M)
    spatial_neg = np.exp(-1j * 2.0 * np.pi * (-mu_cand) * m_idx / M)

    A_pos = np.abs(np.mean(q_analytic * spatial_pos, axis=1))
    A_neg = np.abs(np.mean(q_analytic * spatial_neg, axis=1))

    if np.mean(A_pos) > np.mean(A_neg):
        mu_star = int(mu_cand)
        A_envelope = A_pos
    else:
        mu_star = int(M - mu_cand) if mu_cand > 0 else 0
        A_envelope = A_neg

    # -------------------------------------------------------------------------
    # Step 6: Exponential Fitting & Intrinsic Growth Rate Extraction
    # -------------------------------------------------------------------------
    t_window = np.arange(start_idx, N, dtype=np.float64) * float(downsample_interval)

    if fit_range is not None:
        i0_fit, i1_fit = fit_range
        t_fit = t_window[i0_fit:i1_fit]
        env_fit = A_envelope[i0_fit:i1_fit]
    else:
        # Exclude edge transient artifacts of finite-window Hilbert transform (5% margins)
        margin = max(1, int(0.05 * N_win))
        if N_win > 2 * margin + 2:
            t_fit = t_window[margin : N_win - margin]
            env_fit = A_envelope[margin : N_win - margin]
        else:
            t_fit = t_window
            env_fit = A_envelope

    valid = np.isfinite(env_fit) & (env_fit > 0)
    if np.count_nonzero(valid) < 2:
        growth_fit = 0.0
    else:
        coeffs = np.polyfit(t_fit[valid], np.log(env_fit[valid]), 1)
        growth_fit = float(coeffs[0])

    net_growth_rate = float(growth_fit + (1.0 / tau_rad))

    return mu_star, net_growth_rate, growth_fit, A_envelope


# =============================================================================
# Appendix D: LLRF Open-Loop Transfer Function, Bode & Nyquist Stability Analysis
# =============================================================================
def analyze_llrf_nyquist(
    f_rf: float,
    Q_L: float,
    tau_d: float,
    G_dfb: float = 0.25,
    phi_dfb_deg: float = 0.0,
    detune_hz: float = 0.0,
    m_avg: int = 1,
    T_b: Optional[float] = None,
    K_p: float = 0.0,
    K_i: float = 0.0,
    m_pi: int = 80,
    tau_pi: Optional[float] = None,
    f_min_hz: float = 1e2,
    f_max_hz: float = 5e6,
    num_points: int = 2000,
) -> dict:
    """
    Appendix D LLRF open-loop transfer function L(j*omega_m), Bode response,
    and Nyquist stability margin calculation.

    Open-loop transfer function in rotating frame around omega_rf:
        H_cav(omega_m) = (1 - j*tan(psi)) / (1 + j*(omega_m / omega_half - tan(psi)))
        L(j*omega_m)   = [C_dfb(omega_m) + C_pi(omega_m)] * H_cav(omega_m)

    where:
        omega_half = pi * f_rf / Q_L  (cavity half-bandwidth [rad/s])
        tan(psi)   = 2 * Q_L * detune_hz / f_rf
        C_dfb      = G_dfb * H_avg(omega_m) * exp(-j * (omega_m * tau_d + phi_dfb))

    Returns
    -------
    result : dict
        Dictionary containing frequency grid `f_m_hz`, complex loop responses
        `L_pos` (omega_m > 0), `L_neg` (omega_m < 0), `L_unit` (G_dfb = 1),
        phase crossover frequency `f_pi_hz`, maximum stable gain `G_dfb_max`,
        gain margin `GM_dB`, gain crossover frequency `f_gc_hz`, phase margin `PM_deg`,
        and boolean `is_stable`.
    """
    w_rf = 2.0 * np.pi * f_rf
    w_half = w_rf / (2.0 * Q_L)
    f_half_hz = w_half / (2.0 * np.pi)
    tan_psi = 2.0 * Q_L * detune_hz / f_rf

    f_m_hz = np.logspace(np.log10(f_min_hz), np.log10(f_max_hz), num_points)
    w_m = 2.0 * np.pi * f_m_hz

    # Baseband cavity transfer function (normalized to 1 at DC omega_m = 0)
    H_cav_pos = (1.0 - 1j * tan_psi) / (1.0 + 1j * (w_m / w_half - tan_psi))
    H_cav_neg = (1.0 - 1j * tan_psi) / (1.0 + 1j * (-w_m / w_half - tan_psi))

    # Moving average sampler filter
    if T_b is not None and m_avg > 1:
        x_avg = 0.5 * w_m * m_avg * T_b
        H_avg_pos = (np.sin(x_avg) / np.maximum(1e-30, x_avg)) * np.exp(
            -1j * 0.5 * w_m * (m_avg - 1) * T_b
        )
    else:
        H_avg_pos = np.ones_like(w_m, dtype=np.complex128)
    H_avg_neg = np.conj(H_avg_pos)

    phi_rad = np.deg2rad(phi_dfb_deg)
    H_delay_pos = H_avg_pos * np.exp(-1j * (w_m * tau_d + phi_rad))
    H_delay_neg = H_avg_neg * np.exp(+1j * (w_m * tau_d - phi_rad))

    # Unit-gain DFB open-loop response (for ultimate gain G_dfb_max)
    L_unit_pos = H_delay_pos * H_cav_pos
    L_unit_neg = H_delay_neg * H_cav_neg

    # Optional PI branch
    if (K_p != 0.0 or K_i != 0.0) and T_b is not None:
        Ts_pi = m_pi * T_b
        t_d_pi = tau_pi if tau_pi is not None else tau_d
        x_pi = 0.5 * w_m * m_pi * T_b
        H_avg_pi = (np.sin(x_pi) / np.maximum(1e-30, x_pi)) * np.exp(
            -1j * 0.5 * w_m * (m_pi - 1) * T_b
        )
        C_pi_pos = (
            (K_p + K_i / (1j * w_m * Ts_pi))
            * H_avg_pi
            * np.exp(-1j * w_m * t_d_pi)
        )
        C_pi_neg = (
            (K_p + K_i / (-1j * w_m * Ts_pi))
            * np.conj(H_avg_pi)
            * np.exp(+1j * w_m * t_d_pi)
        )
    else:
        C_pi_pos = 0.0
        C_pi_neg = 0.0

    L_pos = (G_dfb * H_delay_pos + C_pi_pos) * H_cav_pos
    L_neg = (G_dfb * H_delay_neg + C_pi_neg) * H_cav_neg

    # Phase crossover frequency f_pi (where angle(L) = -pi) and Max Stable Gain G_dfb_max
    pha_unit = np.unwrap(np.angle(L_unit_pos))
    idx_u = np.where(pha_unit <= -np.pi)[0]
    if idx_u.size > 0 and idx_u[0] > 0:
        i0 = int(idx_u[0])
        frac = (-np.pi - pha_unit[i0 - 1]) / (pha_unit[i0] - pha_unit[i0 - 1])
        f_pi_unit_hz = float(f_m_hz[i0 - 1] + frac * (f_m_hz[i0] - f_m_hz[i0 - 1]))
        mag_u = float(
            np.abs(L_unit_pos[i0 - 1])
            + frac * (np.abs(L_unit_pos[i0]) - np.abs(L_unit_pos[i0 - 1]))
        )
        G_dfb_max = 1.0 / max(1e-12, mag_u)
    else:
        f_pi_unit_hz = float("nan")
        G_dfb_max = float("inf")

    # Active loop Gain Margin (GM) & Phase Margin (PM)
    pha_act = np.unwrap(np.angle(L_pos))
    mag_act = np.abs(L_pos)
    idx_pc = np.where(pha_act <= -np.pi)[0]
    if idx_pc.size > 0 and idx_pc[0] > 0:
        i0 = int(idx_pc[0])
        frac = (-np.pi - pha_act[i0 - 1]) / (pha_act[i0] - pha_act[i0 - 1])
        f_pi_hz = float(f_m_hz[i0 - 1] + frac * (f_m_hz[i0] - f_m_hz[i0 - 1]))
        mag_pc = float(mag_act[i0 - 1] + frac * (mag_act[i0] - mag_act[i0 - 1]))
        GM_dB = float(-20.0 * np.log10(max(1e-12, mag_pc)))
    else:
        f_pi_hz = f_pi_unit_hz
        mag_pc = 0.0
        GM_dB = float("inf")

    idx_gc = np.where((mag_act[:-1] >= 1.0) & (mag_act[1:] < 1.0))[0]
    if idx_gc.size > 0:
        ig = int(idx_gc[0])
        frac = (1.0 - mag_act[ig]) / (mag_act[ig + 1] - mag_act[ig])
        f_gc_hz = float(f_m_hz[ig] + frac * (f_m_hz[ig + 1] - f_m_hz[ig]))
        pha_gc = float(pha_act[ig] + frac * (pha_act[ig + 1] - pha_act[ig]))
        PM_deg = float(180.0 + np.rad2deg(pha_gc))
    else:
        f_gc_hz = float("nan")
        PM_deg = float("inf")

    is_stable = (GM_dB > 0.0) and (np.isnan(PM_deg) or PM_deg > 0.0)

    return {
        "f_m_hz": f_m_hz,
        "f_half_hz": f_half_hz,
        "L_pos": L_pos,
        "L_neg": L_neg,
        "L_unit_pos": L_unit_pos,
        "L_unit_neg": L_unit_neg,
        "f_pi_hz": f_pi_hz,
        "mag_at_pi": mag_pc,
        "G_dfb_max": G_dfb_max,
        "GM_dB": GM_dB,
        "f_gc_hz": f_gc_hz,
        "PM_deg": PM_deg,
        "is_stable": is_stable,
    }


def plot_llrf_bode_nyquist(
    res: dict,
    gains_sweep: Tuple[float, ...] = (0.2, 0.3, 0.5, 1.0),
    save_path: Optional[str] = None,
) -> None:
    """
    Plot Bode (Magnitude & Phase) and Nyquist diagrams for LLRF stability analysis
    (Appendix D).
    """
    import matplotlib.pyplot as plt

    f_kHz = res["f_m_hz"] * 1e-3
    f_half_kHz = res["f_half_hz"] * 1e-3
    f_pi_kHz = res["f_pi_hz"] * 1e-3
    g_max = res["G_dfb_max"]
    L_unit_pos = res["L_unit_pos"]
    L_unit_neg = res["L_unit_neg"]

    gains = list(gains_sweep)
    if np.isfinite(g_max) and g_max <= 6.0:
        gains.append(round(g_max, 2))
    gains = sorted(set(gains))

    fig = plt.figure(figsize=(13.5, 6.2), dpi=160)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.05, 1.1], hspace=0.28, wspace=0.24)
    ax_mag = fig.add_subplot(gs[0, 0])
    ax_pha = fig.add_subplot(gs[1, 0], sharex=ax_mag)
    ax_nyq = fig.add_subplot(gs[:, 1])

    colors = ["#1f77b4", "#2ca02c", "#ff7f0e", "#9467bd", "#d62728"]

    for idx, g in enumerate(gains):
        col = colors[idx % len(colors)]
        is_crit = np.isfinite(g_max) and np.isclose(g, g_max, atol=0.02)
        lw = 2.2 if is_crit or np.isclose(g, 0.3, atol=0.05) else 1.6
        lbl = f"$G_{{dfb}} = {g:.2f}$ (Crit. $G_{{max}}$)" if is_crit else f"$G_{{dfb}} = {g:.2f}$"

        L_p = g * L_unit_pos
        L_n = g * L_unit_neg
        mag_db = 20.0 * np.log10(np.maximum(1e-12, np.abs(L_p)))
        pha_deg = np.rad2deg(np.unwrap(np.angle(L_p)))

        ax_mag.semilogx(f_kHz, mag_db, color=col, lw=lw, label=lbl)
        ax_pha.semilogx(f_kHz, pha_deg, color=col, lw=lw, label=lbl)

        ax_nyq.plot(np.real(L_p), np.imag(L_p), color=col, lw=lw, label=lbl)
        ax_nyq.plot(np.real(L_n), np.imag(L_n), color=col, lw=1.0, ls=":")

    # Bode formatting
    ax_mag.axhline(0.0, color="k", ls="--", lw=1.1)
    ax_mag.axvline(f_half_kHz, color="gray", ls=":", lw=1.2, label=f"$f_{{1/2}} = {f_half_kHz:.1f}$ kHz")
    if np.isfinite(f_pi_kHz):
        ax_mag.axvline(f_pi_kHz, color="#d62728", ls=":", lw=1.2, label=f"$f_\\pi = {f_pi_kHz:.1f}$ kHz")
    ax_mag.set_ylabel("Magnitude $|L(j\\omega_m)|$ [dB]", fontsize=10.5)
    ax_mag.set_title("LLRF Open-Loop Bode Diagram (Appendix D)", fontsize=11.5, fontweight="bold")
    ax_mag.set_ylim(-42, 22)
    ax_mag.grid(True, which="both", alpha=0.35)
    ax_mag.legend(loc="lower left", fontsize=8.5, ncol=2)

    ax_pha.axhline(-180.0, color="#d62728", ls="--", lw=1.3, label="$-180^\\circ$ Stability Limit")
    if np.isfinite(f_pi_kHz):
        ax_pha.plot(f_pi_kHz, -180.0, "ro", ms=6)
    ax_pha.set_xlabel("Modulation Frequency $f_m = (\\omega - \\omega_{rf}) / 2\\pi$ [kHz]", fontsize=10.5)
    ax_pha.set_ylabel("Phase $\\angle L(j\\omega_m)$ [deg]", fontsize=10.5)
    ax_pha.set_ylim(-340, 20)
    ax_pha.set_xlim(f_kHz[0], f_kHz[-1])
    ax_pha.grid(True, which="both", alpha=0.35)
    ax_pha.legend(loc="lower left", fontsize=8.5)

    # Nyquist formatting
    theta = np.linspace(0, 2.0 * np.pi, 300)
    ax_nyq.plot(np.cos(theta), np.sin(theta), color="gray", ls="--", lw=1.2, label="Unit Circle ($|L|=1$)")
    ax_nyq.axhline(0.0, color="gray", ls=":", lw=0.9)
    ax_nyq.axvline(0.0, color="gray", ls=":", lw=0.9)
    ax_nyq.plot(-1.0, 0.0, "r+", ms=12, mew=2.4, label="Critical Point $(-1, j0)$")
    ax_nyq.plot(-1.0, 0.0, "ro", ms=6, fillstyle="none", mew=1.4)
    ax_nyq.set_aspect("equal", adjustable="box")
    ax_nyq.set_xlim(-2.1, 2.4)
    ax_nyq.set_ylim(-2.1, 2.1)
    ax_nyq.set_xlabel("$\\mathrm{Re}[L(j\\omega_m)]$", fontsize=10.5)
    ax_nyq.set_ylabel("$\\mathrm{Im}[L(j\\omega_m)]$", fontsize=10.5)
    ax_nyq.set_title(
        f"Nyquist Diagram ($f_\\pi={f_pi_kHz:.1f}$ kHz, Max Stable Gain $G_{{max}}={g_max:.2f}$)",
        fontsize=11.5,
        fontweight="bold",
    )
    ax_nyq.grid(True, alpha=0.35)
    ax_nyq.legend(loc="upper right", fontsize=8.5)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=160, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# Unit Tests
# =============================================================================
def _run_unit_tests() -> None:
    """
    Synthetic verification of Appendix C instability mode & growth rate extraction:
      Test 1: Positive coupled-bunch mode (+mu_true = 17)
      Test 2: Negative coupled-bunch mode (-mu_true -> M - mu_true = 416 - 23 = 393)
    """
    rng = np.random.default_rng(42)

    N_samples = 3000
    M_bunches = 416
    downsample_interval = 10.0  # 10 turns per sample
    tau_rad = 11500.0           # radiation damping time in turns
    omega_s = 2.0 * np.pi * 0.0045  # synchrotron tune Q_s = 0.0045

    t_turns = np.arange(N_samples, dtype=np.float64)[:, np.newaxis] * downsample_interval
    m_idx = np.arange(M_bunches, dtype=np.float64)[np.newaxis, :]

    # --- Test 1: Positive mode mu = +17 ---
    mu_pos_true = 17
    growth_fit_true = 1.2e-4  # 1/turn observed growth rate
    net_gr_true = growth_fit_true + 1.0 / tau_rad

    # Physical coupled-bunch wave: q[n, m] = A0 * exp(g * t) * cos(omega_s * t - 2*pi*mu*m/M) + DC
    dc_offset = 0.05
    q_pos = (
        dc_offset
        + 1e-4 * np.exp(growth_fit_true * t_turns)
        * np.cos(omega_s * t_turns - 2.0 * np.pi * mu_pos_true * m_idx / M_bunches)
        + 1e-8 * rng.standard_normal((N_samples, M_bunches))
    )

    mu_star_1, net_gr_1, gr_fit_1, env_1 = extract_tracking_instability(
        q_pos, tau_rad=tau_rad, downsample_interval=downsample_interval
    )

    assert mu_star_1 == mu_pos_true, f"Expected mu_star={mu_pos_true}, got {mu_star_1}"
    assert np.isclose(gr_fit_1, growth_fit_true, rtol=1e-2), (
        f"Expected growth_fit={growth_fit_true:.6e}, got {gr_fit_1:.6e}"
    )
    assert np.isclose(net_gr_1, net_gr_true, rtol=1e-2), (
        f"Expected net_growth_rate={net_gr_true:.6e}, got {net_gr_1:.6e}"
    )
    assert env_1.shape == (N_samples - int(2 * N_samples // 3),), f"Unexpected envelope shape {env_1.shape}"

    # --- Test 2: Negative mode mu = -23 (i.e. M - 23 = 393) ---
    mu_neg_mag = 23
    mu_neg_expected = M_bunches - mu_neg_mag  # 393
    growth_fit_true_2 = 1.8e-4
    net_gr_true_2 = growth_fit_true_2 + 1.0 / tau_rad

    q_neg = (
        -0.02
        + 1e-4 * np.exp(growth_fit_true_2 * t_turns)
        * np.cos(omega_s * t_turns - 2.0 * np.pi * (-mu_neg_mag) * m_idx / M_bunches)
        + 1e-8 * rng.standard_normal((N_samples, M_bunches))
    )

    mu_star_2, net_gr_2, gr_fit_2, env_2 = extract_tracking_instability(
        q_neg, tau_rad=tau_rad, downsample_interval=downsample_interval
    )

    assert mu_star_2 == mu_neg_expected, f"Expected mu_star={mu_neg_expected}, got {mu_star_2}"
    assert np.isclose(gr_fit_2, growth_fit_true_2, rtol=1e-2), (
        f"Expected growth_fit={growth_fit_true_2:.6e}, got {gr_fit_2:.6e}"
    )
    assert np.isclose(net_gr_2, net_gr_true_2, rtol=1e-2), (
        f"Expected net_growth_rate={net_gr_true_2:.6e}, got {net_gr_2:.6e}"
    )

    print("=== All Appendix C Unit Tests Passed ===")
    print(
        f"[Test 1: +mu] True mu={mu_pos_true}, Extracted mu*={mu_star_1} | "
        f"True gr_fit={growth_fit_true:.6e}, Extracted={gr_fit_1:.6e} | "
        f"Net gr={net_gr_1:.6e}"
    )
    print(
        f"[Test 2: -mu] True mu={mu_neg_expected} (-{mu_neg_mag}), Extracted mu*={mu_star_2} | "
        f"True gr_fit={growth_fit_true_2:.6e}, Extracted={gr_fit_2:.6e} | "
        f"Net gr={net_gr_2:.6e}"
    )

    # --- Test 3: Appendix D LLRF Nyquist & Bode Analysis ---
    f_rf = 352.2e6
    Q_L = 12000.0
    tau_d = 1.5e-6  # 1.5 us loop delay
    res_llrf = analyze_llrf_nyquist(
        f_rf=f_rf, Q_L=Q_L, tau_d=tau_d, G_dfb=0.30, phi_dfb_deg=0.0, detune_hz=-20e3
    )
    assert res_llrf["is_stable"] is True, "Expected G_dfb=0.30 to be stable"
    assert res_llrf["G_dfb_max"] > 1.0, f"Unexpected G_dfb_max={res_llrf['G_dfb_max']}"
    print("=== All Appendix D Unit Tests Passed ===")
    print(
        f"[Test 3: LLRF Nyquist] f_half={res_llrf['f_half_hz']*1e-3:.2f} kHz | "
        f"f_pi={res_llrf['f_pi_hz']*1e-3:.2f} kHz | "
        f"Max Stable Gain G_max={res_llrf['G_dfb_max']:.2f} | "
        f"GM={res_llrf['GM_dB']:.2f} dB | Stable={res_llrf['is_stable']}"
    )


if __name__ == "__main__":
    _run_unit_tests()
