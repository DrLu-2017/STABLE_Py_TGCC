"""
Single-Point Multi-Bunch Longitudinal GPU Tracking Engine for STABLE_Py_TGCC.
Faithfully ports MATLAB STABLE_App.runSimulation & MPMBTrackBunchLengthingwoHighQ.m
with live streaming TelemetryFrame callbacks for remote GUI plotting.
"""

from __future__ import annotations
from typing import Callable, Optional, Dict, Any, List
from collections import deque
import ctypes
import time
import math
import numpy as np

from stable_tgcc.common.models import MachineConfig, TelemetryFrame
from stable_tgcc.physics.backend import ArrayBackend
from stable_tgcc.physics.machine_init import (
    truncated_gaussian,
    generate_fill_pattern,
    detune_hc_calc,
    compute_optimal_mc_detune,
    align_hom_to_worst_sideband,
    build_half_machine,
)
from stable_tgcc.physics.cavity_feedback import (
    dfb_set,
    mc_cavity_control,
    voltage_load_calc,
)
from stable_tgcc.physics.analysis import classify_tracking_result


def run_single_tracking(
    cfg: MachineConfig,
    job_id: str = "job_single",
    device_id: int = 0,
    prefer_cpu: bool = False,
    on_telemetry: Optional[Callable[[TelemetryFrame], None]] = None,
    stop_check: Optional[Callable[[], bool]] = None,
    telemetry_updates: int = 40,
) -> TelemetryFrame:
    """
    Execute full multi-bunch longitudinal tracking on GPU (`device_id`) or CPU,
    emitting live `TelemetryFrame` updates via `on_telemetry`.
    """
    cfg.sync_derived_cavity_params()
    backend = ArrayBackend(device_id=device_id, prefer_cpu=prefer_cpu)

    cspeed = 299792458.0
    C = float(cfg.C)
    E0 = float(cfg.E0)
    I0 = float(cfg.I0)
    U0 = float(cfg.U0)
    h = int(cfg.h)
    alpha_c = float(cfg.alpha_c)
    tau_s = float(cfg.tau_s)
    tau_z = float(cfg.tau_z)
    sigma_e0 = float(cfg.sigma_e0)

    V_mc = float(cfg.V_mc)
    Q_mc = float(cfg.Q_mc)
    R_mc = float(cfg.R_mc)
    n_hc = int(cfg.n_hc)
    Q_hc = float(cfg.Q_hc)
    R_hc = float(cfg.R_hc)

    Par_num = max(10, int(cfg.Par_num))
    Track_num = max(10, int(cfg.Track_num))
    Recor_step = max(1, int(cfg.Recor_step))

    # Natural bunch length sigma_t0
    Qs_est = math.sqrt(max(1e-12, alpha_c / (2.0 * math.pi) * h * math.sqrt(max(0.0, (V_mc / E0) ** 2 - (U0 / E0) ** 2))))
    T0_est = C / cspeed
    sigma_t0 = alpha_c * sigma_e0 / (2.0 * math.pi * Qs_est) * T0_est

    # Fill pattern
    pattern = generate_fill_pattern(h, cfg.pattern_mode, cfg.custom_pattern_expr)
    Bun_num = int(np.sum(pattern))
    fillrate = Bun_num / float(h)
    filled_mask = pattern == 1

    # HHC detuning
    if cfg.hc_auto_detune:
        fre_shift_hc = detune_hc_calc(
            I0, n_hc, C, h, U0, V_mc, R_hc, Q_hc, E0, alpha_c, sigma_e0, cfg.use_form_factor_iter, 1.0
        )
    else:
        fre_shift_hc = float(cfg.fre_shift_hc)

    # MC detuning
    if cfg.mc_auto_detune:
        fre_shift_mc = compute_optimal_mc_detune(
            C, h, I0, V_mc, U0, R_mc, Q_mc, R_hc, Q_hc, fre_shift_hc, n_hc
        )
    else:
        fre_shift_mc = float(cfg.fre_shift_mc)

    # Build HALF machine dictionary
    HALF = build_half_machine(
        C, I0, U0, E0, tau_s, tau_z, sigma_t0, sigma_e0, alpha_c,
        h, V_mc, n_hc, R_hc, Q_hc, fillrate, fre_shift_hc,
        Q_mc, R_mc, fre_shift_mc, cfg.use_form_factor_iter,
    )
    HALF["ShortRange_on"] = 1 if cfg.short_range_wake else 0

    # HOMs setup
    has_homs = len(cfg.homs) > 0
    if has_homs:
        Q_hom = np.array([m.Q for m in cfg.homs], dtype=np.float64)
        R_hom = np.array([m.R for m in cfg.homs], dtype=np.float64)
        fre_hom = np.array([m.freq_Hz for m in cfg.homs], dtype=np.float64)
        if cfg.hom_auto_align:
            fre_hom = align_hom_to_worst_sideband(fre_hom, C, h, E0, U0, V_mc, alpha_c)
        wrf_hom = 2.0 * math.pi * fre_hom
        angle_hom = wrf_hom * HALF["Tb"]
        im_factor_hom = 1.0 / np.sqrt(4.0 * (Q_hom ** 2) - 1.0)
        rot_coef_hom = np.sqrt(1.0 - 1.0 / ((2.0 * Q_hom) ** 2))
        n_hom = len(Q_hom)
    else:
        n_hom = 0

    # Feedback setup (PI + DFB)
    PI: Dict[str, Any] = {
        "m": max(1, int(cfg.pi_m)),
        "d": max(1, int(cfg.pi_d)),
        "piIndex": 0,
        "Integral": 0j,
        "RL": HALF["R_mc"],
        "Ts": max(1, int(cfg.pi_m)) * HALF["Tb"],
        "KP": float(cfg.pi_kp),
        "KI": float(cfg.pi_ki),
        "Ig0": HALF["Ig_mc_0"],
        "DIg": 0j,
        "Ig": HALF["Ig_mc_0"],
        "Vg0": HALF["Vg_mc_0"],
        "RoverQ": HALF["R_mc"] / HALF["Q_mc"],
        "Vrf_mc_track": np.zeros(max(1, int(cfg.pi_m)), dtype=np.complex128),
        "Ig_FBid": deque(),
        "enable": bool(cfg.pi_enable),
    }

    if cfg.dfb_enable:
        dfb_gain_val = cfg.get_dfb_gain_for_current(I0)
        dfb_phase_rad = math.radians(cfg.dfb_phase_shift_deg)
        DFB: Optional[Dict[str, Any]] = dfb_set(
            HALF,
            gain=dfb_gain_val,
            phase_shift=dfb_phase_rad,
            delay=cfg.dfb_delay,
            sample_num=cfg.dfb_m,
            every=cfg.dfb_every,
        )
        DFB["enable"] = True
        PI["Ig0"] = DFB["Ig_main0"]
        PI["Ig"] = DFB["Ig_main0"]
        PI["Vg0"] = DFB["Vg0_main"]
    else:
        DFB = None

    # Initialize macro-particles on GPU/CPU
    q1 = truncated_gaussian(1.0, (-3.0, 3.0), (Par_num, 1))
    p1 = truncated_gaussian(1.0, (-3.0, 3.0), (Par_num, 1))
    Q = backend.to_device(np.repeat(q1, Bun_num, axis=1), dtype=np.float32)
    P = backend.to_device(np.repeat(p1, Bun_num, axis=1), dtype=np.float32)

    # Charge & cavity coupling coefficients
    charge = np.ones(h, dtype=np.float64) * pattern
    charge = charge / max(1.0, float(np.sum(charge))) * Bun_num
    qc_arr = charge * pattern * HALF["qc"] / float(Par_num)

    Vb_hc_filled = (qc_arr * HALF["wr_hc"] * R_hc / Q_hc)[filled_mask]
    Vb_mc_filled = (qc_arr * HALF["wr_mc"] * R_mc / Q_mc)[filled_mask]

    V_hc_load_0 = complex(HALF["V_hc_load_0"])
    V_mc_load_0 = complex(HALF["V_mc_load_0"])
    Vg_mc_track_0 = complex(HALF["Vg_mc_init"])

    Dq = 0.05
    rot_decay_hc = 1j * HALF["rot_coef_hc"] - 1.0 / (2.0 * Q_hc)
    TbAng_coef_hc = np.exp(rot_decay_hc * HALF["angle_hc"])
    exp_ang_coef_hc = -rot_decay_hc * HALF["wr_hc"] * sigma_t0

    rot_decay_mc = 1j * HALF["rot_coef_mc"] - 1.0 / (2.0 * Q_mc)
    TbAng_coef_mc = np.exp(rot_decay_mc * HALF["angle_mc"])
    exp_ang_coef_mc = -rot_decay_mc * HALF["wr_mc"] * sigma_t0

    if has_homs:
        Vb_hom_filled = (qc_arr[None, :] * ((wrf_hom * R_hom / Q_hom)[:, None]))[:, filled_mask]
        rot_decay_hom = 1j * rot_coef_hom - 1.0 / (2.0 * Q_hom)
        TbAng_coef_hom = np.exp(rot_decay_hom * angle_hom)
        exp_ang_coef_hom = -rot_decay_hom * wrf_hom * sigma_t0
        V_hom_load_0 = np.zeros(n_hom, dtype=np.complex128)

    wake_kick_coef_filled = (qc_arr * HALF["kick_coef"])[filled_mask]

    # Intra-bunch short-range wake table
    tau_q = np.arange(0.0, 200.0 + 0.5 * Dq, Dq) * sigma_t0
    Wake_inter_hc = (
        -HALF["wr_hc"] * (R_hc / Q_hc)
        * np.exp(-tau_q * HALF["wr_hc"] / (2.0 * Q_hc))
        * (
            np.cos(tau_q * HALF["wr_hc"] * HALF["rot_coef_hc"])
            - HALF["VbImagFactor_hc"] * np.sin(tau_q * HALF["wr_hc"] * HALF["rot_coef_hc"])
        )
    )
    Wake_inter_hc[0] *= 0.5

    Wake_inter_mc = (
        -HALF["wr_mc"] * (R_mc / Q_mc)
        * np.exp(-tau_q * HALF["wr_mc"] / (2.0 * Q_mc))
        * (
            np.cos(tau_q * HALF["wr_mc"] * HALF["rot_coef_mc"])
            - HALF["VbImagFactor_mc"] * np.sin(tau_q * HALF["wr_mc"] * HALF["rot_coef_mc"])
        )
    )
    Wake_inter_mc[0] *= 0.5

    Wake_inter = Wake_inter_hc + Wake_inter_mc
    if has_homs:
        for jj in range(n_hom):
            w_hom = (
                -wrf_hom[jj] * (R_hom[jj] / Q_hom[jj])
                * np.exp(-tau_q * wrf_hom[jj] / (2.0 * Q_hom[jj]))
                * (
                    np.cos(tau_q * wrf_hom[jj] * rot_coef_hom[jj])
                    - im_factor_hom[jj] * np.sin(tau_q * wrf_hom[jj] * rot_coef_hom[jj])
                )
            )
            w_hom[0] *= 0.5
            Wake_inter = Wake_inter + w_hom

    # Recording buffers
    Recor_num = max(1, Track_num // Recor_step)
    record_Q_mean = np.zeros((Recor_num, Bun_num), dtype=np.float64)
    record_Q_std = np.zeros((Recor_num, Bun_num), dtype=np.float64)
    record_P_mean = np.zeros((Recor_num, Bun_num), dtype=np.float64)
    record_P_std = np.zeros((Recor_num, Bun_num), dtype=np.float64)
    record_th = 0

    turns_axis: List[int] = []
    bl_history_ps: List[float] = []
    jit_history_ps: List[float] = []
    es_history: List[float] = []

    ui_interval = max(25, Track_num // max(1, telemetry_updates))
    drift_coef_f32 = np.float32(HALF["drift_coef"])
    radamp_f32 = np.float32(HALF["radampcoef"])
    quanex_f32 = np.float32(HALF["quanexcoef"])
    ploss_f32 = np.float32(HALF["ploss"])
    kick_coef = float(HALF["kick_coef"])

    V_load_hc_last = np.zeros(h, dtype=np.complex128)
    Vc_mc_last = np.zeros(h, dtype=np.complex128)

    t_start = time.perf_counter()
    beam_lost = False

    if backend.mode == "nvcc_cuda":
        lib = backend.cuda_lib
        q_col = np.ascontiguousarray(np.repeat(q1, Bun_num, axis=1).T, dtype=np.float32)
        p_col = np.ascontiguousarray(np.repeat(p1, Bun_num, axis=1).T, dtype=np.float32)
        w_f32 = np.ascontiguousarray(Wake_inter[:2005], dtype=np.float32)
        wc_f32 = np.ascontiguousarray(wake_kick_coef_filled, dtype=np.float32)

        ctx = lib.stable_gpu_init(
            int(backend.device_id),
            int(Par_num),
            int(Bun_num),
            q_col.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            p_col.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            w_f32.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
            len(w_f32),
            wc_f32.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
        )
        sum_re_buf = np.zeros(Bun_num, dtype=np.float64)
        sum_im_buf = np.zeros(Bun_num, dtype=np.float64)
        qm_buf = np.zeros(Bun_num, dtype=np.float64)
        qs_buf = np.zeros(Bun_num, dtype=np.float64)
        pm_buf = np.zeros(Bun_num, dtype=np.float64)
        ps_buf = np.zeros(Bun_num, dtype=np.float64)
        rng_seed = int(time.time_ns() & 0xFFFFFFFFFFFFFFFF)

        try:
            for turn in range(1, Track_num + 1):
                if stop_check is not None and stop_check():
                    break

                max_q = float(lib.stable_gpu_drift_and_max(ctx, float(drift_coef_f32)))
                max_tau_abs = max_q * sigma_t0
                if math.isnan(max_tau_abs) or math.isinf(max_tau_abs) or max_tau_abs > 2800e-12:
                    beam_lost = True
                    break

                # HOMs kick
                if has_homs:
                    for jj in range(n_hom):
                        c_hom = exp_ang_coef_hom[jj]
                        lib.stable_gpu_exp_and_sum(
                            ctx,
                            float(np.real(c_hom)),
                            float(np.imag(c_hom)),
                            sum_re_buf.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                            sum_im_buf.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        )
                        sum_hom = sum_re_buf + 1j * sum_im_buf
                        v_load_cpu = sum_hom * Vb_hom_filled[jj]
                        v_load_all, V_hom_load_0[jj] = voltage_load_calc(
                            V_hom_load_0[jj], v_load_cpu, TbAng_coef_hom[jj], pattern
                        )
                        vk = (v_load_all[filled_mask] * kick_coef).astype(np.complex64)
                        vr = np.ascontiguousarray(np.real(vk), dtype=np.float32)
                        vi = np.ascontiguousarray(np.imag(vk), dtype=np.float32)
                        lib.stable_gpu_apply_cavity_kick(
                            ctx,
                            vr.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                            vi.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                            float(im_factor_hom[jj]),
                        )

                # Radiation damping + quantum excitation + short-range wake
                ok_w = lib.stable_gpu_rad_and_sr_wake(
                    ctx,
                    float(1.0 / Dq),
                    float(min(turn / 5000.0, 1.0)),
                    1 if cfg.short_range_wake else 0,
                    float(radamp_f32),
                    float(quanex_f32),
                    float(ploss_f32),
                    rng_seed,
                    turn,
                )
                if ok_w == 0:
                    beam_lost = True
                    break

                # HHC kick
                lib.stable_gpu_exp_and_sum(
                    ctx,
                    float(np.real(exp_ang_coef_hc)),
                    float(np.imag(exp_ang_coef_hc)),
                    sum_re_buf.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    sum_im_buf.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                )
                sum_hc = sum_re_buf + 1j * sum_im_buf
                v_load_hc_cpu = sum_hc * Vb_hc_filled
                V_load_hc_last, V_hc_load_0 = voltage_load_calc(
                    V_hc_load_0, v_load_hc_cpu, TbAng_coef_hc, pattern
                )
                vk_hc = (V_load_hc_last[filled_mask] * kick_coef).astype(np.complex64)
                vr_hc = np.ascontiguousarray(np.real(vk_hc), dtype=np.float32)
                vi_hc = np.ascontiguousarray(np.imag(vk_hc), dtype=np.float32)
                lib.stable_gpu_apply_cavity_kick(
                    ctx,
                    vr_hc.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                    vi_hc.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                    float(HALF["VbImagFactor_hc"]),
                )

                # MC kick + PI & DFB Feedbacks
                lib.stable_gpu_exp_and_sum(
                    ctx,
                    float(np.real(exp_ang_coef_mc)),
                    float(np.imag(exp_ang_coef_mc)),
                    sum_re_buf.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    sum_im_buf.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                )
                sum_mc = sum_re_buf + 1j * sum_im_buf
                v_load_mc_cpu = sum_mc * Vb_mc_filled
                Vc_mc_last, _, Vg_mc_track_0, _, V_mc_load_0, PI, DFB = mc_cavity_control(
                    PI, DFB, HALF["Vrf_ideal"], Vg_mc_track_0, V_mc_load_0, v_load_mc_cpu, TbAng_coef_mc, pattern
                )
                vk_mc = (Vc_mc_last[filled_mask] * kick_coef).astype(np.complex64)
                vr_mc = np.ascontiguousarray(np.real(vk_mc), dtype=np.float32)
                vi_mc = np.ascontiguousarray(np.imag(vk_mc), dtype=np.float32)
                lib.stable_gpu_apply_cavity_kick(
                    ctx,
                    vr_mc.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                    vi_mc.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                    float(HALF["VbImagFactor_mc"]),
                )

                if turn % Recor_step == 0 and record_th < Recor_num:
                    lib.stable_gpu_mean_std(
                        ctx,
                        qm_buf.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        qs_buf.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        pm_buf.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                        ps_buf.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                    )
                    record_Q_mean[record_th, :] = qm_buf
                    record_Q_std[record_th, :] = qs_buf
                    record_P_mean[record_th, :] = pm_buf
                    record_P_std[record_th, :] = ps_buf
                    record_th += 1

                    turns_axis.append(turn)
                    bl_history_ps.append(float(np.mean(qs_buf) * sigma_t0 * 1e12))
                    w_start = max(0, record_th - 20)
                    jit_cur = float(np.mean(np.std(record_Q_mean[w_start:record_th, :], axis=0)) * sigma_t0 * 1e12)
                    jit_history_ps.append(jit_cur)
                    es_history.append(float(np.mean(ps_buf)))

                if on_telemetry is not None and (turn % ui_interval == 0 or turn == Track_num):
                    elapsed = time.perf_counter() - t_start
                    eta = (elapsed / max(1, turn)) * (Track_num - turn)
                    frame = TelemetryFrame(
                        job_id=job_id,
                        event_type="turn_progress",
                        turn=turn,
                        total_turns=Track_num,
                        elapsed_sec=round(elapsed, 2),
                        eta_sec=round(eta, 1),
                        device_name=backend.device_name,
                        beam_lost=beam_lost,
                        hc_amp_kV=np.minimum(np.abs(V_load_hc_last) * 1e-3, 3000.0).round(2).tolist(),
                        hc_phase_deg=np.rad2deg(np.angle(V_load_hc_last)).round(2).tolist(),
                        mc_amp_kV=np.minimum(np.abs(Vc_mc_last) * 1e-3, 3000.0).round(2).tolist(),
                        mc_phase_deg=np.rad2deg(np.angle(Vc_mc_last)).round(2).tolist(),
                        turns_axis=turns_axis[-200:],
                        bunch_length_ps_history=bl_history_ps[-200:],
                        centroid_jitter_ps_history=jit_history_ps[-200:],
                        energy_spread_history=es_history[-200:],
                        xi_eff=float(HALF["xi_eff"]),
                        psi_hc_deg=float(HALF["psi_hc_deg"]),
                    )
                    on_telemetry(frame)
        finally:
            lib.stable_gpu_free(ctx)
    else:
        for turn in range(1, Track_num + 1):
            if stop_check is not None and stop_check():
                break

            # 1. Longitudinal drift
            Q = Q + P * drift_coef_f32

            # 2. Beam loss check
            max_tau_abs = backend.max_abs(Q) * sigma_t0
            if math.isnan(max_tau_abs) or math.isinf(max_tau_abs) or max_tau_abs > 2800e-12:
                beam_lost = True
                break

            # 3. HOMs kick
            if has_homs:
                for jj in range(n_hom):
                    exp_hom, sum_hom = backend.exp_and_col_sum(exp_ang_coef_hom[jj], Q)
                    v_load_cpu = sum_hom * Vb_hom_filled[jj]
                    v_load_all, V_hom_load_0[jj] = voltage_load_calc(
                        V_hom_load_0[jj], v_load_cpu, TbAng_coef_hom[jj], pattern
                    )
                    P = backend.apply_cavity_kick(
                        P, v_load_all[filled_mask] * kick_coef, exp_hom, float(im_factor_hom[jj])
                    )

            # 4. HHC kick
            exp_hc, sum_hc = backend.exp_and_col_sum(exp_ang_coef_hc, Q)
            v_load_hc_cpu = sum_hc * Vb_hc_filled
            V_load_hc_last, V_hc_load_0 = voltage_load_calc(
                V_hc_load_0, v_load_hc_cpu, TbAng_coef_hc, pattern
            )

            # 5. MC kick + PI & DFB Feedbacks
            exp_mc, sum_mc = backend.exp_and_col_sum(exp_ang_coef_mc, Q)
            v_load_mc_cpu = sum_mc * Vb_mc_filled
            Vc_mc_last, _, Vg_mc_track_0, _, V_mc_load_0, PI, DFB = mc_cavity_control(
                PI, DFB, HALF["Vrf_ideal"], Vg_mc_track_0, V_mc_load_0, v_load_mc_cpu, TbAng_coef_mc, pattern
            )

            # 6. Short-range wake kick
            if cfg.short_range_wake:
                ok_wake, wake_kick = backend.short_range_wake_kick(
                    Q, Dq, Wake_inter, wake_kick_coef_filled, min(turn / 5000.0, 1.0)
                )
                if not ok_wake:
                    beam_lost = True
                    break
            else:
                wake_kick = np.float32(0.0)

            # 7. Radiation damping + quantum excitation + cavity kicks
            P = (
                P
                - radamp_f32 * P
                + quanex_f32 * backend.randn((Par_num, Bun_num))
                - ploss_f32
                + wake_kick
            )
            P = backend.apply_cavity_kick(
                P, V_load_hc_last[filled_mask] * kick_coef, exp_hc, float(HALF["VbImagFactor_hc"])
            )
            P = backend.apply_cavity_kick(
                P, Vc_mc_last[filled_mask] * kick_coef, exp_mc, float(HALF["VbImagFactor_mc"])
            )

            # 8. Record statistics
            if turn % Recor_step == 0 and record_th < Recor_num:
                qm, qs, pm, ps = backend.mean_and_std_cols(Q, P)
                record_Q_mean[record_th, :] = qm
                record_Q_std[record_th, :] = qs
                record_P_mean[record_th, :] = pm
                record_P_std[record_th, :] = ps
                record_th += 1

                turns_axis.append(turn)
                bl_history_ps.append(float(np.mean(qs) * sigma_t0 * 1e12))
                w_start = max(0, record_th - 20)
                jit_cur = float(np.mean(np.std(record_Q_mean[w_start:record_th, :], axis=0)) * sigma_t0 * 1e12)
                jit_history_ps.append(jit_cur)
                es_history.append(float(np.mean(ps)))

            # 9. Emit live telemetry frame
            if on_telemetry is not None and (turn % ui_interval == 0 or turn == Track_num):
                elapsed = time.perf_counter() - t_start
                eta = (elapsed / max(1, turn)) * (Track_num - turn)
                frame = TelemetryFrame(
                    job_id=job_id,
                    event_type="turn_progress",
                    turn=turn,
                    total_turns=Track_num,
                    elapsed_sec=round(elapsed, 2),
                    eta_sec=round(eta, 1),
                    device_name=backend.device_name,
                    beam_lost=beam_lost,
                    hc_amp_kV=np.minimum(np.abs(V_load_hc_last) * 1e-3, 3000.0).round(2).tolist(),
                    hc_phase_deg=np.rad2deg(np.angle(V_load_hc_last)).round(2).tolist(),
                    mc_amp_kV=np.minimum(np.abs(Vc_mc_last) * 1e-3, 3000.0).round(2).tolist(),
                    mc_phase_deg=np.rad2deg(np.angle(Vc_mc_last)).round(2).tolist(),
                    turns_axis=turns_axis[-200:],
                    bunch_length_ps_history=bl_history_ps[-200:],
                    centroid_jitter_ps_history=jit_history_ps[-200:],
                    energy_spread_history=es_history[-200:],
                    xi_eff=float(HALF["xi_eff"]),
                    psi_hc_deg=float(HALF["psi_hc_deg"]),
                )
                on_telemetry(frame)

    backend.synchronize()
    elapsed_total = time.perf_counter() - t_start

    valid_rec = max(1, record_th)
    diag = classify_tracking_result(
        record_Q_mean[:valid_rec, :],
        record_Q_std[:valid_rec, :],
        record_P_mean[:valid_rec, :],
        record_P_std[:valid_rec, :],
        sigma_t0=sigma_t0,
        T0=HALF["T0"],
        tau_z=HALF["tau_z"],
        Recor_step=Recor_step,
        beam_lost=beam_lost,
        xi_eff=float(HALF["xi_eff"]),
        I0=I0,
        fre_shift_hc=fre_shift_hc,
        Par_num=Par_num,
        has_homs=has_homs,
    )

    final_frame = TelemetryFrame(
        job_id=job_id,
        event_type="tracking_done",
        turn=turn if beam_lost else Track_num,
        total_turns=Track_num,
        elapsed_sec=round(elapsed_total, 2),
        eta_sec=0.0,
        device_name=backend.device_name,
        beam_lost=beam_lost,
        hc_amp_kV=np.minimum(np.abs(V_load_hc_last) * 1e-3, 3000.0).round(2).tolist(),
        hc_phase_deg=np.rad2deg(np.angle(V_load_hc_last)).round(2).tolist(),
        mc_amp_kV=np.minimum(np.abs(Vc_mc_last) * 1e-3, 3000.0).round(2).tolist(),
        mc_phase_deg=np.rad2deg(np.angle(Vc_mc_last)).round(2).tolist(),
        turns_axis=turns_axis,
        bunch_length_ps_history=bl_history_ps,
        centroid_jitter_ps_history=jit_history_ps,
        energy_spread_history=es_history,
        stability_code=int(diag["stability_code"]),
        status_label=str(diag["status_label"]),
        growth_rate_s1=float(diag["growth_rate_s1"]),
        dominant_freq_Hz=float(diag["dominant_freq_Hz"]),
        dominant_mode_mu=int(diag["dominant_mode_mu"]),
        mode_spectrum_amp_ps=diag["mode_spectrum_amp_ps"],
        xi_eff=float(HALF["xi_eff"]),
        psi_hc_deg=float(HALF["psi_hc_deg"]),
    )
    if on_telemetry is not None:
        on_telemetry(final_frame)
    return final_frame
