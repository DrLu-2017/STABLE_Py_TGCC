"""
Multithreaded HTTP + Real-Time NDJSON Streaming RPC Server for STABLE_Py_TGCC.
Runs on TGCC node (reached via SSH local port forwarding `ssh -L`) or locally.
Requires zero external web frameworks (built on Python's `http.server.ThreadingHTTPServer`).
"""

from __future__ import annotations
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Dict, Any, Set, Tuple
import json
import socket
import threading
import traceback

from stable_tgcc.common.models import MachineConfig, ScanConfig, TelemetryFrame, ScanPointTelemetry
from stable_tgcc.common.protocol import encode_ndjson_event
from stable_tgcc.physics.backend import detect_available_gpus
from stable_tgcc.physics.machine_init import compute_optimal_mc_detune, calc_bunch_form_factor_eq
from stable_tgcc.physics.analysis import analyze_llrf_nyquist, calculate_d_mode_robinson
from stable_tgcc.physics.tracking_engine import run_single_tracking
from stable_tgcc.physics.multi_gpu_scanner import run_2d_stability_scan


class _JobRegistry:
    """Thread-safe registry of running jobs and cancellation flags."""

    def __init__(self):
        self._lock = threading.Lock()
        self._cancelled: Set[str] = set()
        self._active: Set[str] = set()

    def start_job(self, job_id: str) -> None:
        with self._lock:
            self._cancelled.discard(job_id)
            self._active.add(job_id)

    def finish_job(self, job_id: str) -> None:
        with self._lock:
            self._active.discard(job_id)
            self._cancelled.discard(job_id)

    def cancel_job(self, job_id: str) -> None:
        with self._lock:
            if job_id == "*":
                self._cancelled.update(self._active)
            else:
                self._cancelled.add(job_id)

    def is_cancelled(self, job_id: str) -> bool:
        with self._lock:
            return job_id in self._cancelled

    def active_jobs(self) -> list[str]:
        with self._lock:
            return sorted(self._active)


_REGISTRY = _JobRegistry()


class _STABLERequestHandler(BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.1"

    def log_message(self, format: str, *args: Any) -> None:
        # Suppress noisy per-request console spam
        pass

    def _read_json_body(self) -> Dict[str, Any]:
        length = int(self.headers.get("Content-Length", "0"))
        if length <= 0:
            return {}
        raw = self.rfile.read(length)
        return json.loads(raw.decode("utf-8"))

    def _send_json(self, status_code: int, payload: Dict[str, Any]) -> None:
        body = json.dumps(payload).encode("utf-8")
        self.send_response(status_code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Connection", "close")
        self.end_headers()
        self.wfile.write(body)
        self.wfile.flush()

    def do_GET(self) -> None:
        if self.path.startswith("/api/status"):
            gpus = detect_available_gpus()
            self._send_json(
                200,
                {
                    "ok": True,
                    "hostname": socket.gethostname(),
                    "gpus": gpus,
                    "active_jobs": _REGISTRY.active_jobs(),
                },
            )
        else:
            self._send_json(404, {"ok": False, "error": f"Unknown path: {self.path}"})

    def do_POST(self) -> None:
        try:
            body = self._read_json_body()
            if self.path.startswith("/api/stop"):
                job_id = str(body.get("job_id", "*"))
                _REGISTRY.cancel_job(job_id)
                self._send_json(200, {"ok": True, "cancelled": job_id})
                return

            if self.path.startswith("/api/nyquist"):
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
                self._send_json(200, {"ok": True, "nyquist": res})
                return

            if self.path.startswith("/api/dmode"):
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
                self._send_json(200, {"ok": True, "dmode": dmode, "equilibrium": fp_clean})
                return

            if self.path.startswith("/api/track_stream"):
                self._handle_track_stream(body)
                return

            if self.path.startswith("/api/scan_stream"):
                self._handle_scan_stream(body)
                return

            self._send_json(404, {"ok": False, "error": f"Unknown endpoint: {self.path}"})
        except Exception as exc:
            self._send_json(500, {"ok": False, "error": str(exc), "trace": traceback.format_exc()})

    def _handle_track_stream(self, body: Dict[str, Any]) -> None:
        job_id = str(body.get("job_id", "job_track"))
        device_id = int(body.get("device_id", 0))
        cfg = MachineConfig.from_dict(body.get("machine", {}))

        self.send_response(200)
        self.send_header("Content-Type", "application/x-ndjson; charset=utf-8")
        self.send_header("Cache-Control", "no-cache")
        self.send_header("Connection", "close")
        self.end_headers()

        write_lock = threading.Lock()

        def emit_frame(frame: TelemetryFrame) -> None:
            with write_lock:
                try:
                    self.wfile.write(encode_ndjson_event("telemetry", frame.to_dict()))
                    self.wfile.flush()
                except Exception:
                    _REGISTRY.cancel_job(job_id)

        _REGISTRY.start_job(job_id)
        try:
            run_single_tracking(
                cfg=cfg,
                job_id=job_id,
                device_id=device_id,
                on_telemetry=emit_frame,
                stop_check=lambda: _REGISTRY.is_cancelled(job_id),
            )
        except Exception as exc:
            with write_lock:
                try:
                    self.wfile.write(encode_ndjson_event("error", {"job_id": job_id, "error": str(exc)}))
                    self.wfile.flush()
                except Exception:
                    pass
        finally:
            _REGISTRY.finish_job(job_id)

    def _handle_scan_stream(self, body: Dict[str, Any]) -> None:
        job_id = str(body.get("job_id", "job_scan"))
        cfg = MachineConfig.from_dict(body.get("machine", {}))
        scan_cfg = ScanConfig.from_dict(body.get("scan", {}))

        self.send_response(200)
        self.send_header("Content-Type", "application/x-ndjson; charset=utf-8")
        self.send_header("Cache-Control", "no-cache")
        self.send_header("Connection", "close")
        self.end_headers()

        write_lock = threading.Lock()

        def emit_point(pt: ScanPointTelemetry) -> None:
            with write_lock:
                try:
                    self.wfile.write(encode_ndjson_event("scan_point", pt.to_dict()))
                    self.wfile.flush()
                except Exception:
                    _REGISTRY.cancel_job(job_id)

        def emit_prog(prog: Dict[str, Any]) -> None:
            with write_lock:
                try:
                    self.wfile.write(encode_ndjson_event("scan_progress", prog))
                    self.wfile.flush()
                except Exception:
                    _REGISTRY.cancel_job(job_id)

        _REGISTRY.start_job(job_id)
        try:
            pts = run_2d_stability_scan(
                base_cfg=cfg,
                scan_cfg=scan_cfg,
                job_id=job_id,
                on_point_done=emit_point,
                on_turn_progress=emit_prog,
                stop_check=lambda: _REGISTRY.is_cancelled(job_id),
            )
            with write_lock:
                self.wfile.write(
                    encode_ndjson_event(
                        "scan_done",
                        {"job_id": job_id, "total_completed": len(pts)},
                    )
                )
                self.wfile.flush()
        except Exception as exc:
            with write_lock:
                try:
                    self.wfile.write(encode_ndjson_event("error", {"job_id": job_id, "error": str(exc)}))
                    self.wfile.flush()
                except Exception:
                    pass
        finally:
            _REGISTRY.finish_job(job_id)


class STABLERPCDaemon:
    """Wrapper around `ThreadingHTTPServer` for clean lifecycle management."""

    def __init__(self, host: str = "127.0.0.1", port: int = 18765):
        self.server = ThreadingHTTPServer((host, port), _STABLERequestHandler)
        self.host, self.port = self.server.server_address[:2]
        self._thread: threading.Thread | None = None

    def serve_forever(self) -> None:
        self.server.serve_forever()

    def start_background(self) -> Tuple[str, int]:
        self._thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self._thread.start()
        return str(self.host), int(self.port)

    def shutdown(self) -> None:
        self.server.shutdown()
        self.server.server_close()


def start_daemon_in_background(host: str = "127.0.0.1", port: int = 0) -> STABLERPCDaemon:
    """Start an in-process local RPC daemon on an ephemeral or fixed port."""
    daemon = STABLERPCDaemon(host=host, port=port)
    daemon.start_background()
    return daemon
