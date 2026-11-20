"""
TGCC Topaze / Irene A100 Batch Scheduler (`ccc_msub` / SLURM) & Direct SSH Stdio-RPC Runner.
Supports multi-GPU sharded execution (`--shard-index k --num-shards N`) where N independent
A100 GPU jobs execute N grid points simultaneously and stream live results to `telemetry.jsonl`.
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Any, List, Optional, Callable
import argparse
import json
import os
import shutil
import socket
import subprocess
import sys
import tempfile
import threading

from stable_tgcc.common.models import (
    MachineConfig,
    ScanConfig,
    TGCCConnectionProfile,
    TelemetryFrame,
    ScanPointTelemetry,
)
from stable_tgcc.common.protocol import encode_ndjson_event
from stable_tgcc.physics.backend import detect_available_gpus
from stable_tgcc.physics.machine_init import compute_optimal_mc_detune, calc_bunch_form_factor_eq
from stable_tgcc.physics.analysis import analyze_llrf_nyquist, calculate_d_mode_robinson
from stable_tgcc.physics.tracking_engine import run_single_tracking
from stable_tgcc.physics.multi_gpu_scanner import run_2d_stability_scan


def generate_ccc_msub_script(
    profile: TGCCConnectionProfile,
    job_id: str,
    job_dir_remote: str,
    job_type: str = "scan",
    num_gpus: int = 1,
    include_account: bool = True,
    shard_index: int = 0,
    num_shards: int = 1,
    qos_override: Optional[str] = None,
    walltime_override: Optional[int] = None,
) -> str:
    """
    Generate a TGCC Topaze / Irene A100 `ccc_msub` job script (`job_submit_s{shard_index}.sh`).
    When `num_shards > 1`, each shard requests 1 A100 GPU (`#MSUB -c 32` without `#MSUB -x`)
    so Topaze can schedule all shards immediately onto `mix` / `mix-` A100 nodes in parallel!
    Supports `#MSUB -Q <qos>` (`test`, `normal`, `long`) per TGCC priority rules (`ccc_mqinfo`).
    """
    n_gpus = max(1, min(4, int(num_gpus)))
    cores = 32 * n_gpus
    qos_val = (qos_override if qos_override is not None else getattr(profile, "ccc_qos", "normal")).strip()
    walltime = max(300, int(walltime_override if walltime_override is not None else profile.walltime_seconds))
    if qos_val == "test":
        walltime = min(walltime, 1800)  # TGCC 'test' QoS max walltime is 30 min (1800s) for high-priority instant start
    work_root = str(Path(job_dir_remote).parents[1].as_posix())
    acct_line = f"#MSUB -A {profile.ccc_project}\n" if (include_account and profile.ccc_project.strip()) else ""
    qos_line = f"#MSUB -Q {qos_val}\n" if qos_val and qos_val not in ("default", "none") else ""
    excl_line = "#MSUB -x\n" if n_gpus == 4 else ""
    shard_tag = f"_s{shard_index}" if num_shards > 1 else ""

    script = f"""#!/bin/bash
#MSUB -r STABLE_{job_id}{shard_tag}
#MSUB -q {profile.ccc_partition}
{qos_line}{acct_line}#MSUB -n 1
#MSUB -c {cores}
{excl_line}#MSUB -T {walltime}
#MSUB -m scratch,work
#MSUB -o {job_dir_remote}/job{shard_tag}_%I.out
#MSUB -e {job_dir_remote}/job{shard_tag}_%I.err

set -e
source /etc/profile >/dev/null 2>&1 || true
module load cuda >/dev/null 2>&1 || true
for mod in pytorch-gpu pytorch cupy python3 python; do
    module load $mod >/dev/null 2>&1 || true
done
export PYTHONPATH="{work_root}:${{PYTHONPATH}}"
cd "{job_dir_remote}"

echo "[TGCC Topaze A100] Starting STABLE job {job_id} (GPU shard {shard_index}/{num_shards}, {job_type}) on $(hostname) at $(date)"
nvidia-smi --query-gpu=index,name,memory.total --format=csv || true

# Enable NVIDIA Multi-Process Service (MPS) if permitted so K worker processes share the A100's 108 SMs concurrently
export CUDA_MPS_PIPE_DIRECTORY="/tmp/mps_${{SLURM_JOB_ID:-$$}}_s{shard_index}"
export CUDA_MPS_LOG_DIRECTORY="/tmp/mps_log_${{SLURM_JOB_ID:-$$}}_s{shard_index}"
mkdir -p "$CUDA_MPS_PIPE_DIRECTORY" "$CUDA_MPS_LOG_DIRECTORY" 2>/dev/null || true
nvidia-cuda-mps-control -d >/dev/null 2>&1 || true
trap 'echo quit | nvidia-cuda-mps-control >/dev/null 2>&1 || true; rm -rf "$CUDA_MPS_PIPE_DIRECTORY" "$CUDA_MPS_LOG_DIRECTORY" 2>/dev/null || true' EXIT

({profile.python_env_cmd} -m stable_tgcc.server.ccc_scheduler \\
    --job-file "{job_dir_remote}/job_request.json" \\
    --telemetry-file "{job_dir_remote}/telemetry.jsonl" \\
    --shard-index {shard_index} \\
    --num-shards {num_shards}) || \\
(python3 -m stable_tgcc.server.ccc_scheduler \\
    --job-file "{job_dir_remote}/job_request.json" \\
    --telemetry-file "{job_dir_remote}/telemetry.jsonl" \\
    --shard-index {shard_index} \\
    --num-shards {num_shards})
"""
    return script


def _make_atomic_point_claimer(claim_dir: Path, candidate_indices: List[int]) -> Callable[[], Optional[int]]:
    """
    Lock-free, POSIX/Lustre/Windows-atomic sequential work queue using `os.mkdir`.
    Ensures all parallel worker processes claim points in strict row-by-row sequential order
    (0, 1, 2, ...) instead of jumping +K indices ahead when an unstable point finishes early!
    """
    claim_dir.mkdir(parents=True, exist_ok=True)
    cursor = [0]

    def _claim() -> Optional[int]:
        while cursor[0] < len(candidate_indices):
            k = candidate_indices[cursor[0]]
            cursor[0] += 1
            try:
                os.mkdir(str(claim_dir / f"pt_{k}"))
                return k
            except FileExistsError:
                continue
            except OSError:
                continue
        return None

    return _claim


def run_batch_job_cli(
    job_file: str | Path,
    telemetry_file: str | Path,
    shard_index: int = 0,
    num_shards: int = 1,
    is_subworker: bool = False,
) -> None:
    """
    Executed on TGCC A100 compute node by `ccc_msub`:
    Reads `job_request.json` and appends NDJSON lines atomically (`O_APPEND`) into `telemetry.jsonl`.
    When `workers_per_gpu > 1` (e.g. 1 A100 GPU runs 5 points concurrently), spawns `K = workers_per_gpu`
    independent OS worker processes on the same A100 GPU so each point has its own Python GIL and CPU cores!
    """
    req: Dict[str, Any] = json.loads(Path(job_file).read_text(encoding="utf-8"))
    job_id = str(req.get("job_id", "batch_job"))
    job_type = str(req.get("job_type", "single"))
    cfg = MachineConfig.from_dict(req.get("machine", {}))
    t_path = Path(telemetry_file)
    t_path.parent.mkdir(parents=True, exist_ok=True)

    write_lock = threading.Lock()
    fd = os.open(str(t_path), os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o644)
    try:
        def write_event(ev_type: str, data: Dict[str, Any]) -> None:
            line_bytes = encode_ndjson_event(ev_type, data)
            with write_lock:
                os.write(fd, line_bytes)

        if job_type == "scan":
            scan_cfg = ScanConfig.from_dict(req.get("scan", {}))
            k_per_gpu = max(1, min(32, int(scan_cfg.workers_per_gpu)))

            if not is_subworker:
                write_event(
                    "job_started",
                    {
                        "job_id": job_id,
                        "shard_index": shard_index,
                        "num_shards": num_shards,
                        "workers_per_gpu": k_per_gpu,
                        "total_concurrent": num_shards * k_per_gpu,
                        "hostname": socket.gethostname(),
                        "gpus": detect_available_gpus(),
                    },
                )

            if not is_subworker and k_per_gpu > 1:
                # Spawn K independent OS processes sharing this A100 GPU (zero Python GIL contention!)
                # All K processes pull sequentially from a shared atomic claim directory so they never jump +K ahead!
                claim_dir = t_path.parent / f".claims_{job_id}_s{shard_index}"
                shutil.rmtree(claim_dir, ignore_errors=True)
                claim_dir.mkdir(parents=True, exist_ok=True)

                threads_per_proc = str(max(1, 32 // k_per_gpu))
                child_env = os.environ.copy()
                child_env["OMP_NUM_THREADS"] = threads_per_proc
                child_env["MKL_NUM_THREADS"] = threads_per_proc
                child_env["OPENBLAS_NUM_THREADS"] = threads_per_proc
                child_env["STABLE_CLAIM_DIR"] = str(claim_dir)

                procs = []
                for _w_idx in range(k_per_gpu):
                    cmd = [
                        sys.executable,
                        "-m",
                        "stable_tgcc.server.ccc_scheduler",
                        "--job-file",
                        str(job_file),
                        "--telemetry-file",
                        str(telemetry_file),
                        "--shard-index",
                        str(shard_index),
                        "--num-shards",
                        str(num_shards),
                        "--subworker",
                    ]
                    procs.append(subprocess.Popen(cmd, env=child_env))

                for p in procs:
                    p.wait()

                shutil.rmtree(claim_dir, ignore_errors=True)
                i0_grid = scan_cfg.get_i0_grid_mA()
                df_grid = scan_cfg.get_df_grid_kHz()
                tot_pts = len(i0_grid) * len(df_grid)
                shard_pts_count = sum(1 for pt_k in range(tot_pts) if (pt_k % max(1, num_shards)) == shard_index)
                write_event(
                    "shard_done",
                    {
                        "job_id": job_id,
                        "shard_index": shard_index,
                        "num_shards": num_shards,
                        "shard_completed": shard_pts_count,
                    },
                )
                if num_shards <= 1:
                    write_event("scan_done", {"job_id": job_id, "total_completed": shard_pts_count})
                return

            # Single worker process executing its assigned slice on the A100 GPU
            scan_cfg_local = ScanConfig.from_dict(req.get("scan", {}))
            scan_cfg_local.num_gpus = 1
            scan_cfg_local.workers_per_gpu = 1

            claim_dir_env = os.environ.get("STABLE_CLAIM_DIR", "").strip()
            next_fn = None
            if claim_dir_env:
                tot_pts = len(scan_cfg_local.get_i0_grid_mA()) * len(scan_cfg_local.get_df_grid_kHz())
                cand = [k for k in range(tot_pts) if (k % max(1, num_shards)) == shard_index]
                next_fn = _make_atomic_point_claimer(Path(claim_dir_env), cand)

            pts = run_2d_stability_scan(
                base_cfg=cfg,
                scan_cfg=scan_cfg_local,
                job_id=job_id,
                on_point_done=lambda pt: write_event("scan_point", pt.to_dict()),
                on_turn_progress=lambda pr: write_event("scan_progress", pr),
                shard_index=shard_index,
                num_shards=num_shards,
                next_index_fn=next_fn,
            )
            if not is_subworker:
                write_event(
                    "shard_done",
                    {
                        "job_id": job_id,
                        "shard_index": shard_index,
                        "num_shards": num_shards,
                        "shard_completed": len(pts),
                    },
                )
                if num_shards <= 1:
                    write_event("scan_done", {"job_id": job_id, "total_completed": len(pts)})
        else:
            write_event(
                "job_started",
                {
                    "job_id": job_id,
                    "shard_index": shard_index,
                    "num_shards": num_shards,
                    "hostname": socket.gethostname(),
                    "gpus": detect_available_gpus(),
                },
            )
            run_single_tracking(
                cfg=cfg,
                job_id=job_id,
                device_id=0,
                on_telemetry=lambda frame: write_event("telemetry", frame.to_dict()),
            )
    finally:
        os.close(fd)


def run_stdio_rpc(action: str) -> None:
    """
    Execute an RPC action directly over SSH stdin/stdout without TCP port forwarding.
    """
    act = action.strip().lower()
    if act == "status":
        gpus = detect_available_gpus()
        out = {
            "ok": True,
            "hostname": socket.gethostname(),
            "gpus": gpus,
            "active_jobs": [],
        }
        sys.stdout.write(json.dumps(out) + "\n")
        sys.stdout.flush()
        return

    raw_in = sys.stdin.read()
    body: Dict[str, Any] = json.loads(raw_in) if raw_in.strip() else {}

    if act == "nyquist":
        cfg = MachineConfig.from_dict(body.get("machine", {}))
        cfg.sync_derived_cavity_params()
        df_mc = (
            compute_optimal_mc_detune(
                cfg.C, cfg.h, cfg.I0, cfg.V_mc, cfg.U0, cfg.R_mc, cfg.Q_mc, cfg.R_hc, cfg.Q_hc, cfg.fre_shift_hc, cfg.n_hc
            )
            if cfg.mc_auto_detune
            else cfg.fre_shift_mc
        )
        res = analyze_llrf_nyquist(
            C=cfg.C,
            h=cfg.h,
            Q_mc=cfg.Q_mc,
            fre_shift_mc=df_mc,
            pi_enable=cfg.pi_enable,
            pi_kp=cfg.pi_kp,
            pi_ki=cfg.pi_ki,
            pi_m=cfg.pi_m,
            pi_d=cfg.pi_d,
            dfb_enable=cfg.dfb_enable,
            dfb_gain=cfg.get_dfb_gain_for_current(cfg.I0),
            dfb_phase_shift_deg=cfg.dfb_phase_shift_deg,
            dfb_delay=cfg.dfb_delay,
            dfb_m=cfg.dfb_m,
            include_detune=bool(body.get("include_detune", True)),
        )
        sys.stdout.write(json.dumps({"ok": True, "nyquist": res}) + "\n")
        sys.stdout.flush()
        return

    if act == "dmode":
        cfg = MachineConfig.from_dict(body.get("machine", {}))
        cfg.sync_derived_cavity_params()
        fp = calc_bunch_form_factor_eq(
            "eval", cfg.I0, cfg.V_mc, cfg.U0, cfg.E0, cfg.C, cfg.h,
            cfg.alpha_c, cfg.sigma_e0, cfg.n_hc, cfg.R_hc, cfg.Q_hc,
            cfg.fre_shift_hc, cfg.use_form_factor_iter,
        )
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
            F_n=float(fp["F_hc"]),
        )
        fp_clean = {
            k: ({"real": v.real, "imag": v.imag} if isinstance(v, complex) else v)
            for k, v in fp.items()
        }
        sys.stdout.write(json.dumps({"ok": True, "dmode": dmode, "equilibrium": fp_clean}) + "\n")
        sys.stdout.flush()
        return

    write_lock = threading.Lock()

    def emit_line(ev_type: str, payload: Dict[str, Any]) -> None:
        with write_lock:
            sys.stdout.buffer.write(encode_ndjson_event(ev_type, payload))
            sys.stdout.buffer.flush()

    if act == "track_stream":
        job_id = str(body.get("job_id", "job_track"))
        device_id = int(body.get("device_id", 0))
        cfg = MachineConfig.from_dict(body.get("machine", {}))
        run_single_tracking(
            cfg=cfg,
            job_id=job_id,
            device_id=device_id,
            on_telemetry=lambda fr: emit_line("telemetry", fr.to_dict()),
        )
        return

    if act == "scan_stream":
        job_id = str(body.get("job_id", "job_scan"))
        cfg = MachineConfig.from_dict(body.get("machine", {}))
        scan_cfg = ScanConfig.from_dict(body.get("scan", {}))
        shard_idx = int(body.get("shard_index", 0))
        num_sh = int(body.get("num_shards", 1))
        is_sub = bool(body.get("is_subworker", False))

        k_workers = max(1, min(16, int(scan_cfg.workers_per_gpu)))
        if not is_sub:
            emit_line(
                "job_started",
                {
                    "job_id": job_id,
                    "shard_index": 0,
                    "num_shards": 1,
                    "workers_per_gpu": k_workers,
                    "total_concurrent": k_workers,
                    "hostname": socket.gethostname(),
                    "gpus": detect_available_gpus(),
                },
            )

        if not is_sub and k_workers > 1:
            claim_dir = Path(tempfile.mkdtemp(prefix=f"stable_claims_{job_id}_"))
            threads_per_proc = str(max(1, 32 // k_workers))
            child_env = os.environ.copy()
            child_env["OMP_NUM_THREADS"] = threads_per_proc
            child_env["MKL_NUM_THREADS"] = threads_per_proc
            child_env["OPENBLAS_NUM_THREADS"] = threads_per_proc
            child_env["STABLE_CLAIM_DIR"] = str(claim_dir)

            procs = []
            for _w_idx in range(k_workers):
                sub_payload = json.dumps({
                    "job_id": job_id,
                    "machine": body.get("machine", {}),
                    "scan": body.get("scan", {}),
                    "shard_index": 0,
                    "num_shards": 1,
                    "is_subworker": True,
                }).encode("utf-8")
                p = subprocess.Popen(
                    [sys.executable, "-m", "stable_tgcc.server.ccc_scheduler", "--stdio-rpc", "scan_stream"],
                    stdin=subprocess.PIPE,
                    stdout=sys.stdout,
                    stderr=sys.stderr,
                    env=child_env,
                )
                if p.stdin is not None:
                    p.stdin.write(sub_payload)
                    p.stdin.close()
                procs.append(p)

            for p in procs:
                p.wait()

            shutil.rmtree(claim_dir, ignore_errors=True)
            tot_pts = len(scan_cfg.get_i0_grid_mA()) * len(scan_cfg.get_df_grid_kHz())
            emit_line("scan_done", {"job_id": job_id, "total_completed": tot_pts})
            return

        scan_cfg_local = ScanConfig.from_dict(body.get("scan", {}))
        scan_cfg_local.num_gpus = 1
        scan_cfg_local.workers_per_gpu = 1
        claim_dir_env = os.environ.get("STABLE_CLAIM_DIR", "").strip()
        next_fn = None
        if claim_dir_env:
            tot_pts = len(scan_cfg_local.get_i0_grid_mA()) * len(scan_cfg_local.get_df_grid_kHz())
            cand = [k for k in range(tot_pts) if (k % max(1, num_sh)) == shard_idx]
            next_fn = _make_atomic_point_claimer(Path(claim_dir_env), cand)

        pts = run_2d_stability_scan(
            base_cfg=cfg,
            scan_cfg=scan_cfg_local,
            job_id=job_id,
            on_point_done=lambda pt: emit_line("scan_point", pt.to_dict()),
            on_turn_progress=lambda pr: emit_line("scan_progress", pr),
            shard_index=shard_idx,
            num_shards=num_sh,
            next_index_fn=next_fn,
        )
        if not is_sub:
            emit_line("scan_done", {"job_id": job_id, "total_completed": len(pts)})
        return


def main() -> None:
    parser = argparse.ArgumentParser(description="TGCC A100 Compute Node Batch & Stdio-RPC Runner")
    parser.add_argument("--job-file", help="Path to job_request.json")
    parser.add_argument("--telemetry-file", help="Path to output telemetry.jsonl")
    parser.add_argument("--shard-index", type=int, default=0, help="Shard index (0 .. num_shards-1)")
    parser.add_argument("--num-shards", type=int, default=1, help="Total number of parallel GPU shards")
    parser.add_argument("--subworker", action="store_true", help="Internal flag for per-GPU multi-process subworker")
    parser.add_argument("--stdio-rpc", help="Execute stdio RPC action (status|nyquist|dmode|track_stream|scan_stream)")
    args = parser.parse_args()
    if args.stdio_rpc:
        run_stdio_rpc(args.stdio_rpc)
    elif args.job_file and args.telemetry_file:
        run_batch_job_cli(args.job_file, args.telemetry_file, args.shard_index, args.num_shards, args.subworker)
    else:
        parser.error("Specify either --stdio-rpc or both --job-file and --telemetry-file")


if __name__ == "__main__":
    main()
