"""
Multi-A100 GPU Parallel 2D Stability Map Scanner for STABLE_Py_TGCC.
Distributes (I0, Df_HC) grid points across all allocated NVIDIA A100 GPUs in parallel,
streaming BOTH live per-turn progress AND immediate per-point completion results back to the GUI.
"""

from __future__ import annotations
from typing import Callable, Optional, List, Tuple, Dict, Any
from concurrent.futures import ThreadPoolExecutor, as_completed
import copy
import time

from stable_tgcc.common.models import MachineConfig, ScanConfig, ScanPointTelemetry, TelemetryFrame
from stable_tgcc.physics.backend import detect_available_gpus
from stable_tgcc.physics.machine_init import calc_bunch_form_factor_eq, compute_optimal_mc_detune
from stable_tgcc.physics.analysis import calculate_d_mode_robinson
from stable_tgcc.physics.tracking_engine import run_single_tracking


def _evaluate_single_grid_point(
    base_cfg: MachineConfig,
    scan_cfg: ScanConfig,
    job_id: str,
    pt_idx: int,
    total_pts: int,
    i0_mA: float,
    df_kHz: float,
    gpu_id: int,
    on_turn_progress: Optional[Callable[[Dict[str, Any]], None]] = None,
    stop_check: Optional[Callable[[], bool]] = None,
) -> ScanPointTelemetry:
    """Worker routine executing a single (I0, Df_HC) point on assigned A100 `gpu_id`."""
    t0 = time.perf_counter()
    cfg = copy.deepcopy(base_cfg)
    cfg.I0 = i0_mA * 1e-3
    cfg.fre_shift_hc = df_kHz * 1e3
    cfg.hc_auto_detune = False
    cfg.sync_derived_cavity_params()

    fp_eq = calc_bunch_form_factor_eq(
        "eval",
        cfg.I0,
        cfg.V_mc,
        cfg.U0,
        cfg.E0,
        cfg.C,
        cfg.h,
        cfg.alpha_c,
        cfg.sigma_e0,
        cfg.n_hc,
        cfg.R_hc,
        cfg.Q_hc,
        cfg.fre_shift_hc,
        cfg.use_form_factor_iter,
    )
    psi_hc_deg = float(fp_eq["psi_hc_deg"])
    xi_eff = float(fp_eq["xi_eff"])
    f_hc_eq = float(fp_eq["F_hc"])

    # Emit immediate point-start notification so GUI shows active blinking marker at (Df_hc, I0) right away
    if on_turn_progress is not None:
        on_turn_progress({
            "job_id": job_id,
            "point_index": pt_idx,
            "total_points": total_pts,
            "I0_mA": i0_mA,
            "Df_hc_kHz": df_kHz,
            "psi_hc_deg": round(psi_hc_deg, 2),
            "xi_eff": round(xi_eff, 4),
            "turn": 0,
            "total_turns": cfg.Track_num,
            "gpu_id": gpu_id,
            "device_name": f"GPU:{gpu_id}",
        })

    # 1. Millisecond D-mode Robinson analytical pre-screening filter (Superconducting Q_hc >= 1e5 only)
    # Normal-conducting cavities (e.g. SOLEIL II Q_hc ~ 3.1e4) always use full particle tracking.
    if scan_cfg.enable_dmode_prescreen and cfg.Q_hc >= 1.0e5:
        df_mc = (
            compute_optimal_mc_detune(
                cfg.C, cfg.h, cfg.I0, cfg.V_mc, cfg.U0, cfg.R_mc, cfg.Q_mc, cfg.R_hc, cfg.Q_hc, cfg.fre_shift_hc, cfg.n_hc
            )
            if cfg.mc_auto_detune
            else cfg.fre_shift_mc
        )
        dmode = calculate_d_mode_robinson(
            I0=cfg.I0,
            E_GeV=cfg.E0 * 1e-9,
            C_m=cfg.C,
            alpha_c=cfg.alpha_c,
            tau_z_s=max(cfg.tau_s, cfg.tau_z),
            U0_keV=cfg.U0 * 1e-3,
            Vrf_MV=cfg.V_mc * 1e-6,
            R1_Q1=cfg.R_mc / max(1.0, cfg.Q_mc),
            Q1=cfg.Q_mc,
            df_r1_Hz=df_mc,
            n_harmonic=cfg.n_hc,
            Rn_Qn=cfg.R_hc / max(1.0, cfg.Q_hc),
            Qn=cfg.Q_hc,
            df_rn_Hz=cfg.fre_shift_hc,
            h=cfg.h,
            F_n=f_hc_eq,
        )
        if dmode["applicable"] and dmode["is_unstable"] and dmode["growth_rate_s1"] > 2.0:
            return ScanPointTelemetry(
                job_id=job_id,
                point_index=pt_idx,
                total_points=total_pts,
                I0_mA=i0_mA,
                Df_hc_kHz=df_kHz,
                psi_hc_deg=round(psi_hc_deg, 2),
                xi_eff=round(xi_eff, 4),
                stability_code=4,
                status_label="D-mode Robinson (Pre-screened)",
                growth_rate_s1=round(float(dmode["growth_rate_s1"]), 2),
                dominant_freq_Hz=round(float(dmode["f_Dmode_Hz"]), 1),
                dominant_mode_mu=0,
                mean_bunch_length_ps=round(float(fp_eq["sigma_t_eq"]) * 1e12, 2),
                centroid_jitter_ps=99.9,
                prescreened_dmode=True,
                gpu_id=gpu_id,
                elapsed_sec=round(time.perf_counter() - t0, 3),
            )

    def _handle_turn_frame(fr: TelemetryFrame) -> None:
        if on_turn_progress is not None and fr.event_type == "turn_progress":
            on_turn_progress({
                "job_id": job_id,
                "point_index": pt_idx,
                "total_points": total_pts,
                "I0_mA": i0_mA,
                "Df_hc_kHz": df_kHz,
                "psi_hc_deg": round(fr.psi_hc_deg, 2),
                "xi_eff": round(fr.xi_eff, 4),
                "turn": fr.turn,
                "total_turns": fr.total_turns,
                "elapsed_sec": fr.elapsed_sec,
                "eta_sec": fr.eta_sec,
                "gpu_id": gpu_id,
                "device_name": fr.device_name,
                "hc_amp_kV": fr.hc_amp_kV,
                "hc_phase_deg": fr.hc_phase_deg,
                "mc_amp_kV": fr.mc_amp_kV,
                "mc_phase_deg": fr.mc_phase_deg,
            })

    # 2. Full GPU macro-particle tracking on assigned A100 GPU
    res = run_single_tracking(
        cfg=cfg,
        job_id=f"{job_id}_p{pt_idx}",
        device_id=gpu_id,
        on_telemetry=_handle_turn_frame,
        stop_check=stop_check,
        telemetry_updates=20,
    )
    mean_bl = res.bunch_length_ps_history[-1] if res.bunch_length_ps_history else (fp_eq["sigma_t_eq"] * 1e12)
    mean_jit = res.centroid_jitter_ps_history[-1] if res.centroid_jitter_ps_history else 0.0

    return ScanPointTelemetry(
        job_id=job_id,
        point_index=pt_idx,
        total_points=total_pts,
        I0_mA=i0_mA,
        Df_hc_kHz=df_kHz,
        psi_hc_deg=round(res.psi_hc_deg, 2),
        xi_eff=round(res.xi_eff, 4),
        stability_code=res.stability_code,
        status_label=res.status_label,
        growth_rate_s1=round(res.growth_rate_s1, 2),
        dominant_freq_Hz=round(res.dominant_freq_Hz, 1),
        dominant_mode_mu=res.dominant_mode_mu,
        mean_bunch_length_ps=round(mean_bl, 2),
        centroid_jitter_ps=round(mean_jit, 2),
        prescreened_dmode=False,
        gpu_id=gpu_id,
        elapsed_sec=round(time.perf_counter() - t0, 2),
    )


def run_2d_stability_scan(
    base_cfg: MachineConfig,
    scan_cfg: ScanConfig,
    job_id: str = "scan_2d",
    on_point_done: Optional[Callable[[ScanPointTelemetry], None]] = None,
    on_turn_progress: Optional[Callable[[Dict[str, Any]], None]] = None,
    stop_check: Optional[Callable[[], bool]] = None,
    shard_index: int = 0,
    num_shards: int = 1,
    next_index_fn: Optional[Callable[[], Optional[int]]] = None,
) -> List[ScanPointTelemetry]:
    """
    Execute a 2D (I0 x Df_HC) stability map scan across multiple A100 GPUs.
    Streams both live per-turn progress (`on_turn_progress`) and per-point completion (`on_point_done`).
    If `next_index_fn` is provided, dynamically pulls the next sequential point index (0..total_pts-1)
    so fast workers never jump far ahead across rows.
    """
    i0_grid = scan_cfg.get_i0_grid_mA()
    df_grid = scan_cfg.get_df_grid_kHz()

    all_points: List[Tuple[int, float, float]] = []
    idx = 0
    for i0 in i0_grid:
        for df in df_grid:
            idx += 1
            all_points.append((idx, i0, df))

    total_pts = len(all_points)
    results: List[ScanPointTelemetry] = []

    if next_index_fn is not None:
        while True:
            if stop_check is not None and stop_check():
                break
            k = next_index_fn()
            if k is None or k < 0 or k >= total_pts:
                break
            pt_idx, i0_mA, df_kHz = all_points[k]
            res = _evaluate_single_grid_point(
                base_cfg,
                scan_cfg,
                job_id,
                pt_idx,
                total_pts,
                i0_mA,
                df_kHz,
                0,
                on_turn_progress,
                stop_check,
            )
            results.append(res)
            if on_point_done is not None:
                on_point_done(res)
        results.sort(key=lambda r: r.point_index)
        return results

    shard_points = [pt for k, pt in enumerate(all_points) if (k % max(1, num_shards)) == shard_index]

    gpus = detect_available_gpus()
    is_cpu_only = (len(gpus) == 1 and gpus[0].get("backend") == "numpy")
    if is_cpu_only:
        # On CPU, run 1 point at a time so Point #1 finishes as fast as possible!
        active_gpus = 1
        max_workers = 1
    else:
        avail_gpu_count = max(1, len(gpus))
        active_gpus = min(max(1, int(scan_cfg.num_gpus)), avail_gpu_count)
        max_workers = max(1, active_gpus * max(1, int(scan_cfg.workers_per_gpu)))

    results: List[ScanPointTelemetry] = []

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = {}
        for k, (pt_idx, i0_mA, df_kHz) in enumerate(shard_points):
            if stop_check is not None and stop_check():
                break
            gpu_id = k % active_gpus
            fut = pool.submit(
                _evaluate_single_grid_point,
                base_cfg,
                scan_cfg,
                job_id,
                pt_idx,
                total_pts,
                i0_mA,
                df_kHz,
                gpu_id,
                on_turn_progress,
                stop_check,
            )
            futures[fut] = pt_idx

        for fut in as_completed(futures):
            res = fut.result()
            results.append(res)
            if on_point_done is not None:
                on_point_done(res)

    results.sort(key=lambda r: r.point_index)
    return results
