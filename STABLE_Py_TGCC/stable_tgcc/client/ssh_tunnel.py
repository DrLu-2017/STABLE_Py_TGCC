"""
SSH Stdio-RPC, SFTP Auto-Deployment & `ccc_msub` Batch Dispatcher for STABLE_Py_TGCC.
"""

from __future__ import annotations
from pathlib import Path
from typing import Callable, Optional, Dict, Any, Tuple, List
import hashlib
import http.client
import io
import json
import re
import tarfile
import threading
import time
import urllib.request

try:
    import paramiko
    _HAS_PARAMIKO = True
except Exception:
    paramiko = None
    _HAS_PARAMIKO = False

from stable_tgcc.common.models import (
    MachineConfig,
    ScanConfig,
    TGCCConnectionProfile,
    TelemetryFrame,
    ScanPointTelemetry,
)
from stable_tgcc.common.protocol import decode_ndjson_line
from stable_tgcc.server.rpc_daemon import start_daemon_in_background, STABLERPCDaemon
from stable_tgcc.server.ccc_scheduler import generate_ccc_msub_script


class STABLEClientConnector:
    """
    Manages the connection between the Client GUI and either:
      - Local embedded RPC server (`profile.mode == "local"`), or
      - Remote TGCC Topaze/Irene A100 cluster over native SSH stdio (`profile.mode == "ssh_tgcc"`).
    """

    def __init__(self):
        self._local_daemon: Optional[STABLERPCDaemon] = None
        self._ssh_client: Any = None
        self._endpoint_host: str = "127.0.0.1"
        self._endpoint_port: int = 0
        self._remote_work_dir: str = "/ccc/work/cont003/soleil/zhaolu/STABLE_Py_TGCC"
        self._remote_python_cmd: str = "python3"
        self._active_profile: Optional[TGCCConnectionProfile] = None
        self._last_batch_id: str = ""
        self._last_batch_ids: List[str] = []
        self._stop_requested: bool = False
        self.active_mode: str = "disconnected"

    def connect(
        self,
        profile: TGCCConnectionProfile,
        ssh_secret: str = "",
        otp_code: str = "",
        log_fn: Optional[Callable[[str], None]] = None,
        step_fn: Optional[Callable[[str, int], None]] = None,
    ) -> Dict[str, Any]:
        """Connect to Local Server or Remote TGCC SSH Cluster and return `/api/status` info."""
        self.disconnect()
        self._active_profile = profile
        self._stop_requested = False

        if profile.mode == "local":
            if step_fn:
                step_fn("Starting local embedded physics server...", 35)
            self._local_daemon = start_daemon_in_background("127.0.0.1", 0)
            self._endpoint_host = "127.0.0.1"
            self._endpoint_port = self._local_daemon.port
            self.active_mode = "local"
            if log_fn:
                log_fn(f"[Local Server] Started embedded RPC daemon on 127.0.0.1:{self._endpoint_port}")
            if step_fn:
                step_fn("Local physics server ready.", 100)
            return self.get_status()

        if not _HAS_PARAMIKO or paramiko is None:
            raise RuntimeError("paramiko is required for SSH connections to TGCC.")

        if step_fn:
            step_fn(f"[1/4] Authenticating SSH to {profile.ssh_username}@{profile.ssh_host}:{profile.ssh_port}...", 18)
        if log_fn:
            log_fn(f"[SSH] Connecting to {profile.ssh_username}@{profile.ssh_host}:{profile.ssh_port}...")

        client = paramiko.SSHClient()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())

        if profile.auth_method == "private_key" and profile.private_key_path:
            client.connect(
                hostname=profile.ssh_host,
                port=int(profile.ssh_port),
                username=profile.ssh_username,
                key_filename=profile.private_key_path,
                passphrase=ssh_secret if ssh_secret else None,
                timeout=20.0,
                banner_timeout=30.0,
            )
        else:
            pwd_combined = f"{ssh_secret}{otp_code}" if otp_code else ssh_secret
            client.connect(
                hostname=profile.ssh_host,
                port=int(profile.ssh_port),
                username=profile.ssh_username,
                password=pwd_combined,
                look_for_keys=True,
                allow_agent=True,
                timeout=20.0,
                banner_timeout=30.0,
            )

        transport = client.get_transport()
        if transport is not None:
            transport.set_keepalive(25)

        self._ssh_client = client
        self.active_mode = "ssh_tgcc"

        if step_fn:
            step_fn("[2/4] Verifying remote work directory & Python environment...", 45)
        raw_dir = profile.remote_work_dir.strip() or "/ccc/work/cont003/soleil/zhaolu/STABLE_Py_TGCC"
        user_py = profile.python_env_cmd.strip() or "python3"
        init_out = self._exec_remote(
            'source /etc/profile >/dev/null 2>&1 || true; '
            f'WD=$(eval echo "{raw_dir}"); '
            'mkdir -p "$WD"; '
            'echo "WD=$WD"; '
            'echo "SHA=$(cat "$WD/.deploy_sha256" 2>/dev/null || true)"; '
            f'({user_py} -c "import numpy, scipy; print(\'PY_OK\')" 2>/dev/null || true)'
        )
        work_dir = f"/ccc/work/cont003/soleil/{profile.ssh_username}/STABLE_Py_TGCC"
        remote_sha = ""
        for ln in init_out.splitlines():
            if ln.startswith("WD=") and ln[3:].strip().startswith("/"):
                work_dir = ln[3:].strip()
            elif ln.startswith("SHA="):
                remote_sha = ln[4:].strip()
        self._remote_work_dir = work_dir

        if "PY_OK" in init_out:
            self._remote_python_cmd = f"source /etc/profile >/dev/null 2>&1 || true; {user_py}"
        else:
            self._remote_python_cmd = (
                "source /etc/profile >/dev/null 2>&1 || true; "
                "(module load python >/dev/null 2>&1 || true); "
                "python3"
            )

        # Synchronize remote Python package (fast SHA256 check + 1-stream tar.gz upload)
        if step_fn:
            step_fn("[3/4] Synchronizing STABLE_Py_TGCC package to Topaze...", 70)
        self.deploy_package_to_tgcc(profile, log_fn=log_fn, force=False, known_remote_sha=remote_sha)

        if log_fn:
            log_fn(f"[SSH Connected] {profile.ssh_username}@{profile.ssh_host} -> {work_dir} (Stdio-RPC ready)")

        if step_fn:
            step_fn("[4/4] Querying remote node status...", 90)
        return self.get_status()

    def disconnect(self) -> None:
        self._stop_requested = True
        if self._ssh_client is not None:
            try:
                self._ssh_client.close()
            except Exception:
                pass
            self._ssh_client = None
        if self._local_daemon is not None:
            try:
                self._local_daemon.shutdown()
            except Exception:
                pass
            self._local_daemon = None
        self.active_mode = "disconnected"

    def _exec_remote(self, cmd: str, timeout: float = 30.0) -> str:
        if self._ssh_client is None:
            raise RuntimeError("Not connected to remote SSH server.")
        _, stdout, stderr = self._ssh_client.exec_command(cmd, timeout=timeout)
        out = stdout.read().decode("utf-8", errors="replace")
        err = stderr.read().decode("utf-8", errors="replace")
        return out if out.strip() else err

    def _run_remote_stdio_rpc_json(self, action: str, payload: Optional[Dict[str, Any]] = None, timeout: float = 30.0) -> Dict[str, Any]:
        """Run a fast JSON RPC call (`status`, `nyquist`, `dmode`) on TGCC over SSH stdin/stdout."""
        if self._ssh_client is None:
            raise RuntimeError("Not connected to TGCC SSH server.")
        cmd = (
            f'cd "{self._remote_work_dir}" && '
            f'export PYTHONPATH="{self._remote_work_dir}:${{PYTHONPATH}}" && '
            f'{self._remote_python_cmd} -m stable_tgcc.server.ccc_scheduler --stdio-rpc {action}'
        )
        stdin, stdout, stderr = self._ssh_client.exec_command(cmd, timeout=timeout)
        if payload is not None:
            stdin.write(json.dumps(payload))
            stdin.flush()
        stdin.channel.shutdown_write()

        out_lines = stdout.read().decode("utf-8", errors="replace").strip().splitlines()
        err_text = stderr.read().decode("utf-8", errors="replace").strip()
        for line in reversed(out_lines):
            line_s = line.strip()
            if line_s.startswith("{") and line_s.endswith("}"):
                return json.loads(line_s)
        raise RuntimeError(f"Remote RPC '{action}' failed: {err_text or out_lines}")

    def deploy_package_to_tgcc(
        self,
        profile: TGCCConnectionProfile,
        log_fn: Optional[Callable[[str], None]] = None,
        force: bool = True,
        known_remote_sha: Optional[str] = None,
    ) -> str:
        """
        Synchronize `stable_tgcc/` and `run_server.py` to `profile.remote_work_dir`.
        Uses SHA-256 change detection + single-stream `.tar.gz` pipe (~0.15s) with SFTP fallback.
        """
        if self._ssh_client is None:
            raise RuntimeError("Connect to TGCC via SSH first before deploying code.")

        local_root = Path(__file__).resolve().parents[2]
        file_list: List[Tuple[Path, str]] = []
        hasher = hashlib.sha256()
        for rel_path in sorted(
            [
                Path("run_server.py"),
                Path("requirements_tgcc_server.txt"),
                *Path(local_root / "stable_tgcc").rglob("*.py"),
            ],
            key=lambda p: p.as_posix(),
        ):
            src = local_root / rel_path if not rel_path.is_absolute() else rel_path
            if not src.exists():
                continue
            rel = src.relative_to(local_root).as_posix()
            raw_bytes = src.read_bytes()
            hasher.update(rel.encode("utf-8"))
            hasher.update(raw_bytes)
            file_list.append((src, rel))

        local_sha = hasher.hexdigest()
        work_dir = self._remote_work_dir
        if not work_dir.startswith("/"):
            raw_dir = profile.remote_work_dir.strip() or "/ccc/work/cont003/soleil/zhaolu/STABLE_Py_TGCC"
            work_dir = self._exec_remote(f'eval echo "{raw_dir}"').strip().splitlines()[-1].strip()
            self._remote_work_dir = work_dir

        if not force:
            rem_sha = known_remote_sha
            if rem_sha is None:
                rem_sha = self._exec_remote(f'cat "{work_dir}/.deploy_sha256" 2>/dev/null || true').strip()
            if rem_sha == local_sha:
                msg = f"Remote TGCC package already up-to-date ({len(file_list)} files, SHA256 verified)"
                if log_fn:
                    log_fn(f"[Fast Sync] {msg}")
                return msg

        # Pack all files into an in-memory .tar.gz (~35 KB) and extract in 1 single SSH stream
        buf = io.BytesIO()
        with tarfile.open(fileobj=buf, mode="w:gz") as tf:
            for src, rel in file_list:
                tf.add(str(src), arcname=rel)
        tar_bytes = buf.getvalue()

        try:
            cmd = f'mkdir -p "{work_dir}" && tar -xzf - -C "{work_dir}" && echo "{local_sha}" > "{work_dir}/.deploy_sha256" && echo "TAR_OK"'
            stdin, stdout, stderr = self._ssh_client.exec_command(cmd, timeout=30.0)
            stdin.write(tar_bytes)
            stdin.flush()
            stdin.channel.shutdown_write()
            out_s = stdout.read().decode("utf-8", errors="replace")
            if "TAR_OK" in out_s:
                msg = f"Deployed {len(file_list)} Python files ({len(tar_bytes)//1024} KB bundle) to {profile.ssh_host}:{work_dir}"
                if log_fn:
                    log_fn(f"[Fast Deploy] {msg}")
                return msg
        except Exception:
            pass

        # Fallback: standard SFTP file-by-file upload
        subdirs = sorted({Path(rel).parent.as_posix() for _, rel in file_list if "/" in rel})
        mkdir_targets = " ".join(f'"{work_dir}/{sd}"' for sd in subdirs)
        self._exec_remote(f'mkdir -p "{work_dir}" {mkdir_targets}')

        sftp = self._ssh_client.open_sftp()
        uploaded = 0
        try:
            for src, rel in file_list:
                sftp.put(str(src), f"{work_dir}/{rel}")
                uploaded += 1
        finally:
            sftp.close()
        self._exec_remote(f'echo "{local_sha}" > "{work_dir}/.deploy_sha256"')

        msg = f"Deployed {uploaded} Python files to {profile.ssh_host}:{work_dir}"
        if log_fn:
            log_fn(f"[SFTP Deploy] {msg}")
        return msg

    def query_tgcc_cluster_info(self) -> str:
        """Query `hostname`, `ccc_mqinfo` (QoS), `ccc_compuse` (Fairshare), `ccc_mpinfo`, `ccc_mpp`, and SLURM `squeue` on Topaze."""
        if self._ssh_client is None:
            status = self.get_status()
            return json.dumps(status, indent=2)
        cmd = (
            'source /etc/profile >/dev/null 2>&1 || true; '
            f'echo "=== TGCC Topaze Node & Work Dir ==="; hostname; ls -ld "{self._remote_work_dir}"; '
            'echo "=== Available QoS & Limits (ccc_mqinfo) ==="; ccc_mqinfo 2>/dev/null || true; '
            'echo "=== Project Fairshare & Allocation (ccc_compuse) ==="; ccc_compuse 2>/dev/null || true; '
            'echo "=== A100 Partition Status (ccc_mpinfo) ==="; ccc_mpinfo -P a100 2>/dev/null || sinfo -p a100 2>/dev/null || true; '
            'echo "=== User Jobs (ccc_mpp) ==="; ccc_mpp -u $USER 2>/dev/null || ccc_mpstat -u $USER 2>/dev/null || true; '
            'echo "=== SLURM Queue Priority, Reason & Est. Start (squeue) ==="; squeue -u $USER -o "%.10i %.8P %.9q %.20j %.8T %.9M %.9l %.10Q %.20S %.25R" 2>/dev/null || true'
        )
        return self._exec_remote(cmd)

    def clean_stale_stable_jobs(self, log_fn: Optional[Callable[[str], None]] = None) -> str:
        """
        Find and cancel any leftover pending/running `STABLE_*` batch jobs belonging to `$USER` on Topaze
        so old queued jobs never block new submissions.
        """
        if self._ssh_client is None:
            return "Not connected to TGCC via SSH."
        cmd = (
            'source /etc/profile >/dev/null 2>&1 || true; '
            'ids=$(squeue -u $USER -h -o "%i %j" 2>/dev/null | awk \'$2 ~ /^STABLE_/ {print $1}\' | tr \'\\n\' \' \'); '
            'if [ -z "$ids" ]; then '
            '  ids=$(ccc_mpp -u $USER 2>/dev/null | awk \'$0 ~ /STABLE_/ {print $3}\' | tr \'\\n\' \' \'); '
            'fi; '
            'if [ -n "$ids" ]; then '
            '  ccc_mdel $ids 2>/dev/null || scancel $ids 2>/dev/null || true; '
            '  echo "CLEANED_IDS=$ids"; '
            'else '
            '  echo "CLEANED_IDS=NONE"; '
            'fi; '
            'ccc_mpp -u $USER 2>/dev/null || true'
        )
        out = self._exec_remote(cmd).strip()
        cleaned = "NONE"
        for ln in out.splitlines():
            if ln.startswith("CLEANED_IDS="):
                cleaned = ln.split("=", 1)[1].strip()
        if cleaned and cleaned != "NONE" and log_fn:
            log_fn(f"[Queue Cleanup] Cancelled previous STABLE batch job(s): {cleaned}")
        return out

    def _ensure_local_connected(self) -> Tuple[str, int]:
        if self.active_mode == "disconnected" or self._endpoint_port <= 0:
            self.connect(TGCCConnectionProfile(mode="local"))
        return self._endpoint_host, self._endpoint_port

    def get_status(self) -> Dict[str, Any]:
        if self.active_mode == "ssh_tgcc" and self._ssh_client is not None:
            return self._run_remote_stdio_rpc_json("status")
        host, port = self._ensure_local_connected()
        url = f"http://{host}:{port}/api/status"
        with urllib.request.urlopen(url, timeout=8.0) as resp:
            return json.loads(resp.read().decode("utf-8"))

    def request_json(self, path: str, payload: Dict[str, Any], timeout: float = 25.0) -> Dict[str, Any]:
        if self.active_mode == "ssh_tgcc" and self._ssh_client is not None:
            action = path.rsplit("/", 1)[-1].strip()
            return self._run_remote_stdio_rpc_json(action, payload=payload, timeout=timeout)
        host, port = self._ensure_local_connected()
        url = f"http://{host}:{port}{path}"
        data = json.dumps(payload).encode("utf-8")
        req = urllib.request.Request(url, data=data, headers={"Content-Type": "application/json"}, method="POST")
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read().decode("utf-8"))

    def stop_job(self, job_id: str = "*") -> Dict[str, Any]:
        self._stop_requested = True
        if self.active_mode == "ssh_tgcc" and self._ssh_client is not None:
            all_ids = list(self._last_batch_ids)
            if self._last_batch_id and self._last_batch_id not in all_ids:
                all_ids.append(self._last_batch_id)
            ids_str = " ".join(all_ids)
            mdel_cmd = (
                f"ccc_mdel {ids_str} 2>/dev/null || scancel {ids_str} 2>/dev/null || true; "
                if ids_str
                else ""
            )
            self._exec_remote(
                "source /etc/profile >/dev/null 2>&1 || true; "
                f"{mdel_cmd}"
                "pkill -f 'stable_tgcc.server.ccc_scheduler' 2>/dev/null || true"
            )
            return {"ok": True, "cancelled": job_id, "batch_ids": all_ids}
        try:
            return self.request_json("/api/stop", {"job_id": job_id}, timeout=5.0)
        except Exception as exc:
            return {"ok": False, "error": str(exc)}

    def _submit_ccc_msub_with_fallback(
        self,
        profile: TGCCConnectionProfile,
        job_id: str,
        job_dir: str,
        job_type: str,
        num_gpus: int,
        log_fn: Optional[Callable[[str], None]],
        shard_index: int = 0,
        num_shards: int = 1,
        print_queue: bool = True,
    ) -> str:
        """
        Submit `job_submit.sh` (or `job_submit_s{shard_index}.sh`) via `ccc_msub` on Topaze.
        If `-A <ccc_project>` is rejected by Topaze, automatically retries without `-A`
        so Topaze uses the user's default project account.
        """
        sh_name = f"job_submit_s{shard_index}.sh" if num_shards > 1 else "job_submit.sh"
        req_qos = getattr(profile, "ccc_qos", "normal").strip() or "normal"
        qos_list = [req_qos, "normal"] if req_qos != "normal" else ["normal"]

        for qos_try in qos_list:
            for include_acct in (True, False):
                msub_sh = generate_ccc_msub_script(
                    profile=profile,
                    job_id=job_id,
                    job_dir_remote=job_dir,
                    job_type=job_type,
                    num_gpus=num_gpus,
                    include_account=include_acct,
                    shard_index=shard_index,
                    num_shards=num_shards,
                    qos_override=qos_try,
                )
                sftp = self._ssh_client.open_sftp()
                try:
                    with sftp.file(f"{job_dir}/{sh_name}", "w") as f_sh:
                        f_sh.write(msub_sh)
                finally:
                    sftp.close()

                cmd = (
                    "source /etc/profile >/dev/null 2>&1 || true; "
                    f'cd "{job_dir}" && (ccc_msub {sh_name} 2>&1 || sbatch {sh_name} 2>&1)'
                )
                sub_out = self._exec_remote(cmd).strip()
                if log_fn:
                    shard_tag = f" [GPU Worker {shard_index + 1}/{num_shards}, QoS={qos_try}]" if num_shards > 1 else f" [QoS={qos_try}]"
                    log_fn(f"[ccc_msub{shard_tag}] {sub_out}")

                m = re.search(r"(\d{5,})", sub_out)
                if m and ("Submitted" in sub_out or "Session" in sub_out or "batch job" in sub_out):
                    batch_id = m.group(1)
                    if print_queue:
                        q_info = self._exec_remote(
                            "source /etc/profile >/dev/null 2>&1 || true; "
                            "ccc_mpp -u $USER 2>/dev/null || ccc_mpstat -u $USER 2>/dev/null || true"
                        ).strip()
                        if log_fn and q_info:
                            log_fn(f"[Queue Status (ccc_mpp)]\n{q_info}")
                    return batch_id

        raise RuntimeError(
            f"ccc_msub submission failed on Topaze. Output: {sub_out}\n"
            "Check 'TGCC Account (-A)' in Tab 6 (see ccc_myproject output in Console)."
        )

    def stream_single_tracking(
        self,
        cfg: MachineConfig,
        job_id: str,
        on_frame: Callable[[TelemetryFrame], None],
        on_error: Optional[Callable[[str], None]] = None,
        device_id: int = 0,
        execution_mode: str = "daemon_stream",
        profile: Optional[TGCCConnectionProfile] = None,
        log_fn: Optional[Callable[[str], None]] = None,
    ) -> None:
        """Stream live `TelemetryFrame` objects over SSH Stdio-RPC, `ccc_msub` A100 batch, or local HTTP NDJSON."""
        self._stop_requested = False
        payload = {"job_id": job_id, "device_id": device_id, "machine": cfg.to_dict()}

        if self.active_mode == "ssh_tgcc" and self._ssh_client is not None:
            if execution_mode == "ccc_msub_batch" and profile is not None:
                self._run_ccc_msub_batch_single(cfg, profile, job_id, on_frame, on_error, log_fn)
                return

            try:
                cmd = (
                    f'cd "{self._remote_work_dir}" && '
                    f'export PYTHONPATH="{self._remote_work_dir}:${{PYTHONPATH}}" && '
                    f'{self._remote_python_cmd} -m stable_tgcc.server.ccc_scheduler --stdio-rpc track_stream'
                )
                stdin, stdout, _ = self._ssh_client.exec_command(cmd)
                stdin.write(json.dumps(payload))
                stdin.flush()
                stdin.channel.shutdown_write()

                for raw_line in iter(stdout.readline, ""):
                    if self._stop_requested:
                        break
                    line_s = raw_line.strip()
                    if not line_s.startswith("{"):
                        continue
                    msg = decode_ndjson_line(line_s)
                    ev = msg.get("event")
                    data = msg.get("data", {})
                    if ev == "telemetry":
                        on_frame(TelemetryFrame.from_dict(data))
                    elif ev == "error" and on_error:
                        on_error(str(data.get("error", "Remote tracking error")))
            except Exception as exc:
                if on_error:
                    on_error(str(exc))
            return

        host, port = self._ensure_local_connected()
        conn = http.client.HTTPConnection(host, port, timeout=3600)
        try:
            body = json.dumps(payload).encode("utf-8")
            conn.request("POST", "/api/track_stream", body=body, headers={"Content-Type": "application/json"})
            resp = conn.getresponse()
            while True:
                if self._stop_requested:
                    break
                line = resp.fp.readline()
                if not line:
                    break
                msg = decode_ndjson_line(line)
                ev = msg.get("event")
                data = msg.get("data", {})
                if ev == "telemetry":
                    on_frame(TelemetryFrame.from_dict(data))
                elif ev == "error" and on_error:
                    on_error(str(data.get("error", "Unknown tracking error")))
        except Exception as exc:
            if on_error:
                on_error(str(exc))
        finally:
            conn.close()

    def stream_2d_scan(
        self,
        cfg: MachineConfig,
        scan_cfg: ScanConfig,
        profile: TGCCConnectionProfile,
        job_id: str,
        on_point: Callable[[ScanPointTelemetry], None],
        on_progress: Optional[Callable[[Dict[str, Any]], None]] = None,
        on_done: Optional[Callable[[int], None]] = None,
        on_error: Optional[Callable[[str], None]] = None,
        log_fn: Optional[Callable[[str], None]] = None,
    ) -> None:
        """
        Execute and stream a 2D Stability Map Scan either via:
          1. TGCC `ccc_msub -q a100` batch job + live SSH file tailing (`execution_mode == "ccc_msub_batch"`),
          2. Direct SSH Stdio-RPC stream (`profile.mode == "ssh_tgcc"` and `execution_mode == "daemon_stream"`), or
          3. Local server stream (`profile.mode == "local"`).
        """
        self._stop_requested = False
        if self.active_mode == "ssh_tgcc" and self._ssh_client is not None:
            if scan_cfg.execution_mode == "ccc_msub_batch":
                self._run_ccc_msub_batch_scan(
                    cfg, scan_cfg, profile, job_id, on_point, on_progress, on_done, on_error, log_fn
                )
                return

            try:
                self.clean_stale_stable_jobs(log_fn=log_fn)
                if log_fn:
                    log_fn(
                        f"[2D Scan — Direct Stream] Starting {scan_cfg.workers_per_gpu} parallel worker processes "
                        f"immediately on {profile.ssh_host} (0s queue wait)..."
                    )
                payload = {
                    "job_id": job_id,
                    "machine": cfg.to_dict(),
                    "scan": scan_cfg.to_dict(),
                }
                cmd = (
                    f'cd "{self._remote_work_dir}" && '
                    f'export PYTHONPATH="{self._remote_work_dir}:${{PYTHONPATH}}" && '
                    f'{self._remote_python_cmd} -m stable_tgcc.server.ccc_scheduler --stdio-rpc scan_stream'
                )
                stdin, stdout, _ = self._ssh_client.exec_command(cmd)
                stdin.write(json.dumps(payload))
                stdin.flush()
                stdin.channel.shutdown_write()

                for raw_line in iter(stdout.readline, ""):
                    if self._stop_requested:
                        break
                    line_s = raw_line.strip()
                    if not line_s.startswith("{"):
                        continue
                    msg = decode_ndjson_line(line_s)
                    ev = msg.get("event")
                    data = msg.get("data", {})
                    if ev == "job_started":
                        if on_progress is not None:
                            on_progress({
                                "event": "job_started",
                                "hostname": data.get("hostname"),
                                "gpus": data.get("gpus"),
                                "shard_index": 0,
                                "num_shards": 1,
                                "workers_per_gpu": data.get("workers_per_gpu", scan_cfg.workers_per_gpu),
                            })
                    elif ev == "scan_progress" and on_progress is not None:
                        on_progress(data)
                    elif ev == "scan_point":
                        on_point(ScanPointTelemetry.from_dict(data))
                    elif ev == "scan_done":
                        if on_done:
                            on_done(int(data.get("total_completed", 0)))
                    elif ev == "error" and on_error:
                        on_error(str(data.get("error", "Remote scan error")))
            except Exception as exc:
                if on_error:
                    on_error(str(exc))
            return

        host, port = self._ensure_local_connected()
        conn = http.client.HTTPConnection(host, port, timeout=86400)
        try:
            body = json.dumps({
                "job_id": job_id,
                "machine": cfg.to_dict(),
                "scan": scan_cfg.to_dict(),
            }).encode("utf-8")
            conn.request("POST", "/api/scan_stream", body=body, headers={"Content-Type": "application/json"})
            resp = conn.getresponse()
            while True:
                if self._stop_requested:
                    break
                line = resp.fp.readline()
                if not line:
                    break
                msg = decode_ndjson_line(line)
                ev = msg.get("event")
                data = msg.get("data", {})
                if ev == "scan_progress" and on_progress is not None:
                    on_progress(data)
                elif ev == "scan_point":
                    on_point(ScanPointTelemetry.from_dict(data))
                elif ev == "scan_done":
                    if on_done:
                        on_done(int(data.get("total_completed", 0)))
                elif ev == "error" and on_error:
                    on_error(str(data.get("error", "Unknown scan error")))
        except Exception as exc:
            if on_error:
                on_error(str(exc))
        finally:
            conn.close()

    def _start_batch_monitor_thread(
        self,
        batch_ids: List[str],
        job_dir: str,
        started_flag: List[bool],
        finished_flag: List[bool],
        tail_chan: Any,
        log_fn: Optional[Callable[[str], None]],
        on_error: Optional[Callable[[str], None]],
        on_progress: Optional[Callable[[Dict[str, Any]], None]] = None,
    ) -> None:
        """Background monitor checking `ccc_mpp`, SLURM `squeue` reason, and `job_<batch_id>.err` while `tail -f` streams."""
        grep_pat = "|".join(re.escape(b) for b in batch_ids if b)
        t_submit = time.perf_counter()

        def _watch():
            while not self._stop_requested and not finished_flag[0]:
                time.sleep(5.0)
                if self._stop_requested or finished_flag[0]:
                    break
                try:
                    q_out = self._exec_remote(
                        "source /etc/profile >/dev/null 2>&1 || true; "
                        f"squeue -u $USER -h -o '%i|%P|%T|%M|%l|%S|%R' 2>/dev/null | grep -E '{grep_pat}' || "
                        f"ccc_mpp -u $USER 2>/dev/null | grep -E '{grep_pat}' || true"
                    ).strip()
                    if q_out:
                        if not started_flag[0]:
                            clean_q = re.sub(r"\x1b\[[0-9;]*m", "", q_out)
                            first_ln = clean_q.splitlines()[0].strip()
                            state_s, reason_s, start_s = "PENDING", "Waiting for A100 slot", "N/A"
                            if "|" in first_ln:
                                parts = [p.strip() for p in first_ln.split("|")]
                                if len(parts) >= 7:
                                    state_s = parts[2]
                                    start_s = parts[5]
                                    reason_s = parts[6]
                            wait_sec = int(round(time.perf_counter() - t_submit))
                            if log_fn:
                                log_fn(f"[Topaze Queue {wait_sec}s] {clean_q}")
                            if on_progress is not None:
                                on_progress({
                                    "event": "queue_status",
                                    "batch_ids": batch_ids,
                                    "state": state_s,
                                    "reason": reason_s,
                                    "start_est": start_s,
                                    "wait_sec": wait_sec,
                                    "raw": clean_q,
                                })
                    else:
                        # No shard jobs remain in ccc_mpp queue; give tail 2s to flush remaining lines
                        time.sleep(2.0)
                        if not finished_flag[0]:
                            err_out = self._exec_remote(
                                f'cat "{job_dir}"/job_*.out "{job_dir}"/job_*.err 2>/dev/null | tail -n 40 || true'
                            ).strip()
                            if err_out and log_fn:
                                log_fn(f"[Batch Output]\n{err_out}")
                            if not started_flag[0] and on_error:
                                on_error(f"Batch job(s) {batch_ids} exited before streaming started. Check console log.")
                            try:
                                tail_chan.close()
                            except Exception:
                                pass
                        break
                except Exception:
                    pass

        threading.Thread(target=_watch, daemon=True).start()

    def _run_ccc_msub_batch_single(
        self,
        cfg: MachineConfig,
        profile: TGCCConnectionProfile,
        job_id: str,
        on_frame: Callable[[TelemetryFrame], None],
        on_error: Optional[Callable[[str], None]],
        log_fn: Optional[Callable[[str], None]],
    ) -> None:
        """Submit a single-point tracking job to Topaze `a100` partition via `ccc_msub` and tail `telemetry.jsonl`."""
        try:
            work_dir = self._remote_work_dir
            job_dir = f"{work_dir}/jobs/{job_id}"
            self._exec_remote(f'mkdir -p "{job_dir}" && : > "{job_dir}/telemetry.jsonl"')

            req_json = json.dumps({
                "job_id": job_id,
                "job_type": "single",
                "machine": cfg.to_dict(),
            }, indent=2)
            sftp = self._ssh_client.open_sftp()
            try:
                with sftp.file(f"{job_dir}/job_request.json", "w") as f_req:
                    f_req.write(req_json)
            finally:
                sftp.close()

            # Clean up any leftover pending STABLE jobs from earlier runs so they don't block the queue
            self.clean_stale_stable_jobs(log_fn=log_fn)

            profile.walltime_seconds = profile.compute_effective_walltime_seconds(cfg, None, is_scan=False)
            if log_fn:
                log_fn(
                    f"[Walltime] Allocated #MSUB -T {profile.walltime_seconds}s "
                    f"({profile.walltime_seconds / 3600.0:.1f}h, auto={profile.auto_walltime})"
                )
            batch_id = self._submit_ccc_msub_with_fallback(profile, job_id, job_dir, "single", 1, log_fn)
            self._last_batch_id = batch_id
            self._last_batch_ids = [batch_id]
            if log_fn:
                log_fn(f"[TGCC Topaze] BatchID {batch_id} submitted to partition '{profile.ccc_partition}'...")

            _, stdout, _ = self._ssh_client.exec_command(f'tail -n +1 -f "{job_dir}/telemetry.jsonl"')
            started_flag = [False]
            finished_flag = [False]
            self._start_batch_monitor_thread(
                [batch_id], job_dir, started_flag, finished_flag, stdout.channel, log_fn, on_error
            )
            try:
                for raw_line in iter(stdout.readline, ""):
                    if self._stop_requested:
                        break
                    line_s = raw_line.strip()
                    if not line_s.startswith("{"):
                        continue
                    msg = decode_ndjson_line(line_s)
                    ev = msg.get("event")
                    data = msg.get("data", {})
                    if ev == "job_started":
                        started_flag[0] = True
                        if log_fn:
                            log_fn(f"[A100 Compute Node Active] Host: {data.get('hostname')} | GPUs: {data.get('gpus')}")
                    elif ev == "telemetry":
                        started_flag[0] = True
                        fr = TelemetryFrame.from_dict(data)
                        on_frame(fr)
                        if fr.event_type == "tracking_done":
                            finished_flag[0] = True
                            break
                    elif ev == "error":
                        finished_flag[0] = True
                        if on_error:
                            on_error(str(data.get("error", "Batch error")))
                        break
            finally:
                finished_flag[0] = True
        except Exception as exc:
            if on_error:
                on_error(str(exc))

    def _run_ccc_msub_batch_scan(
        self,
        cfg: MachineConfig,
        scan_cfg: ScanConfig,
        profile: TGCCConnectionProfile,
        job_id: str,
        on_point: Callable[[ScanPointTelemetry], None],
        on_progress: Optional[Callable[[Dict[str, Any]], None]],
        on_done: Optional[Callable[[int], None]],
        on_error: Optional[Callable[[str], None]],
        log_fn: Optional[Callable[[str], None]],
    ) -> None:
        """
        Submit `num_shards` parallel 1-GPU (`#MSUB -q a100 -c 32`) worker jobs on Topaze
        so each A100 GPU tracks 1 point simultaneously (`1 GPU = 1 Point`), and tail
        the shared `telemetry.jsonl` live over SSH.
        """
        try:
            work_dir = self._remote_work_dir
            i0_grid = scan_cfg.get_i0_grid_mA()
            df_grid = scan_cfg.get_df_grid_kHz()
            total_pts = max(1, int(len(i0_grid) * len(df_grid)))
            num_shards = min(max(1, int(scan_cfg.num_gpus)), total_pts, 16)

            # Check if there is already a RUNNING STABLE_scan_* job on Topaze (e.g. topaze7038);
            # if so, re-attach to its telemetry.jsonl instead of killing the active A100 node!
            running_check = self._exec_remote(
                'source /etc/profile >/dev/null 2>&1 || true; '
                'squeue -u $USER -h -t R -o "%i|%j" 2>/dev/null | grep "STABLE_scan_" || true'
            ).strip()
            attached_existing = False
            batch_ids: List[str] = []
            if running_check:
                exist_job_id = ""
                for ln in running_check.splitlines():
                    if "|" in ln:
                        b_part, j_part = [x.strip() for x in ln.split("|", 1)]
                        batch_ids.append(b_part)
                        # j_part is like STABLE_scan_1795017500 or STABLE_scan_1795017500_s0
                        m_j = re.match(r"^STABLE_(scan_\d+)", j_part)
                        if m_j:
                            exist_job_id = m_j.group(1)
                if exist_job_id and batch_ids:
                    cand_dir = f"{work_dir}/jobs/{exist_job_id}"
                    has_tel = self._exec_remote(f'test -f "{cand_dir}/telemetry.jsonl" && echo "YES" || true').strip()
                    if "YES" in has_tel:
                        job_id = exist_job_id
                        job_dir = cand_dir
                        num_shards = max(1, len(batch_ids))
                        attached_existing = True
                        self._last_batch_id = batch_ids[-1]
                        self._last_batch_ids = list(batch_ids)
                        if log_fn:
                            log_fn(
                                f"[Auto-Attach] Re-attached to already RUNNING A100 job(s) {', '.join(batch_ids)} "
                                f"({job_id}) on Topaze! Replaying completed points & streaming live..."
                            )

            if not attached_existing:
                job_dir = f"{work_dir}/jobs/{job_id}"
                self._exec_remote(f'mkdir -p "{job_dir}" && : > "{job_dir}/telemetry.jsonl"')

                req_json = json.dumps({
                    "job_id": job_id,
                    "job_type": "scan",
                    "machine": cfg.to_dict(),
                    "scan": scan_cfg.to_dict(),
                }, indent=2)

                sftp = self._ssh_client.open_sftp()
                try:
                    with sftp.file(f"{job_dir}/job_request.json", "w") as f_req:
                        f_req.write(req_json)
                finally:
                    sftp.close()

                # Clean up any leftover pending STABLE jobs from earlier runs so they don't block the queue
                self.clean_stale_stable_jobs(log_fn=log_fn)

                k_per_gpu = max(1, int(scan_cfg.workers_per_gpu))
                profile.walltime_seconds = profile.compute_effective_walltime_seconds(cfg, scan_cfg, is_scan=True)
                if log_fn:
                    log_fn(
                        f"[Walltime] Allocated #MSUB -T {profile.walltime_seconds}s "
                        f"({profile.walltime_seconds / 3600.0:.1f}h, auto={profile.auto_walltime}) | "
                        f"Dispatch: {num_shards} A100 GPU(s) × {k_per_gpu} Points/GPU = {num_shards * k_per_gpu} concurrent points."
                    )

                for s_idx in range(num_shards):
                    if self._stop_requested:
                        break
                    b_id = self._submit_ccc_msub_with_fallback(
                        profile=profile,
                        job_id=job_id,
                        job_dir=job_dir,
                        job_type="scan",
                        num_gpus=1,
                        log_fn=log_fn,
                        shard_index=s_idx,
                        num_shards=num_shards,
                        print_queue=(s_idx == num_shards - 1),
                    )
                    batch_ids.append(b_id)

                self._last_batch_id = batch_ids[-1] if batch_ids else ""
                self._last_batch_ids = list(batch_ids)
                if log_fn:
                    log_fn(
                        f"[TGCC Topaze] Submitted {len(batch_ids)} parallel A100 worker job(s) "
                        f"(BatchIDs: {', '.join(batch_ids)}) to partition '{profile.ccc_partition}'. "
                        "Streaming telemetry..."
                    )

            _, stdout, _ = self._ssh_client.exec_command(f'tail -n +1 -f "{job_dir}/telemetry.jsonl"')
            started_flag = [False]
            finished_flag = [False]
            self._start_batch_monitor_thread(
                batch_ids, job_dir, started_flag, finished_flag, stdout.channel, log_fn, on_error, on_progress
            )
            points_received = 0
            shards_done = 0
            seen_points = set()
            try:
                for raw_line in iter(stdout.readline, ""):
                    if self._stop_requested:
                        break
                    line_s = raw_line.strip()
                    if not line_s.startswith("{"):
                        continue
                    msg = decode_ndjson_line(line_s)
                    ev = msg.get("event")
                    data = msg.get("data", {})
                    if ev == "job_started":
                        started_flag[0] = True
                        s_i = int(data.get("shard_index", 0))
                        n_s = int(data.get("num_shards", 1))
                        w_k = int(data.get("workers_per_gpu", max(1, int(scan_cfg.workers_per_gpu))))
                        if log_fn:
                            log_fn(
                                f"[A100 GPU {s_i + 1}/{n_s} Active] "
                                f"Host: {data.get('hostname')} | GPUs: {data.get('gpus')} | "
                                f"Running {w_k} concurrent points on this A100 ({n_s * w_k} total parallel points)"
                            )
                        if on_progress is not None:
                            on_progress({
                                "event": "job_started",
                                "hostname": data.get("hostname"),
                                "gpus": data.get("gpus"),
                                "shard_index": s_i,
                                "num_shards": n_s,
                                "workers_per_gpu": w_k,
                            })
                    elif ev == "scan_progress" and on_progress is not None:
                        started_flag[0] = True
                        on_progress(data)
                    elif ev == "scan_point":
                        started_flag[0] = True
                        pt = ScanPointTelemetry.from_dict(data)
                        pt_key = (int(pt.point_index), round(float(pt.I0_mA), 3), round(float(pt.Df_hc_kHz), 3))
                        if pt_key not in seen_points:
                            seen_points.add(pt_key)
                            points_received += 1
                        on_point(pt)
                        if points_received >= total_pts:
                            finished_flag[0] = True
                            if on_done:
                                on_done(points_received)
                            break
                    elif ev == "shard_done":
                        shards_done += 1
                        if log_fn:
                            log_fn(
                                f"[A100 Worker {int(data.get('shard_index', 0)) + 1}/{num_shards}] "
                                f"Finished shard ({data.get('shard_completed', 0)} pts)."
                            )
                        if shards_done >= num_shards:
                            finished_flag[0] = True
                            if on_done:
                                on_done(points_received)
                            break
                    elif ev == "scan_done":
                        shards_done += 1
                        if shards_done >= num_shards or points_received >= total_pts:
                            finished_flag[0] = True
                            if on_done:
                                on_done(max(points_received, int(data.get("total_completed", 0))))
                            break
                    elif ev == "error":
                        if on_error:
                            on_error(str(data.get("error", "Batch error")))
                        shards_done += 1
                        if shards_done >= num_shards:
                            finished_flag[0] = True
                            break
            finally:
                finished_flag[0] = True
        except Exception as exc:
            if on_error:
                on_error(str(exc))
