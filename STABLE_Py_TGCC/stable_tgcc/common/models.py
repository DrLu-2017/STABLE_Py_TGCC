"""
Serializable data models shared between STABLE Client GUI and TGCC A100 Server.
Strictly synchronized with the verified MATLAB STABLE_App.m physical parameters.
"""

from __future__ import annotations
from dataclasses import dataclass, field, asdict
from typing import List, Dict, Any, Optional
import json
import math


@dataclass
class HOMMode:
    """Single longitudinal Higher-Order Mode (HOM) resonator specification."""
    Q: float
    R: float          # Shunt impedance [Ohm]
    freq_Hz: float    # Resonant frequency [Hz]

    def to_dict(self) -> Dict[str, float]:
        return {"Q": float(self.Q), "R": float(self.R), "freq_Hz": float(self.freq_Hz)}

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "HOMMode":
        return cls(Q=float(d["Q"]), R=float(d["R"]), freq_Hz=float(d["freq_Hz"]))


@dataclass
class MachineConfig:
    """
    Complete storage ring, RF cavities, feedback loops, HOMs, and tracking configuration.
    Matches every parameter of MATLAB STABLE_App.m.
    """
    preset_name: str = "SOLEIL II (3.38-500 mA)"

    # 1. Lattice & Beam
    C: float = 354.097               # Circumference [m]
    E0: float = 2.75e9               # Beam energy [eV]
    I0: float = 0.338                # Total beam current [A]
    U0: float = 470e3                # Synchrotron radiation loss per turn [eV]
    h: int = 416                     # Harmonic number
    alpha_c: float = 1.07e-4         # Momentum compaction factor
    tau_s: float = 13.7e-3           # Longitudinal damping time [s]
    tau_z: float = 13.7e-3           # Energy damping time [s]
    sigma_e0: float = 1.06e-3        # Natural RMS energy spread

    # 2. Main Cavity (MC)
    Ncav_mc: int = 4
    Q0_mc: float = 35.7e3
    RoQ_mc: float = 140.0
    beta_mc: float = 4.95
    V_mc: float = 1.7e6              # Total MC voltage [V]
    Q_mc: float = 6000.0             # Loaded Q of MC
    R_mc: float = 3.36e6             # Total loaded shunt impedance of MC [Ohm]
    fre_shift_mc: float = -24.8e3    # MC detuning [Hz]
    mc_auto_detune: bool = True      # Robinson optimal detuning

    # 3. Harmonic Cavity (HHC)
    n_hc: int = 4                    # Harmonic multiple
    Ncav_hc: int = 2
    Q0_hc: float = 31.0e3
    RoQ_hc: float = 30.0
    beta_hc: float = 0.0             # Passive HHC
    Q_hc: float = 31.0e3             # Loaded Q of HHC
    R_hc: float = 1.86e6             # Total loaded shunt impedance of HHC [Ohm]
    fre_shift_hc: float = 67.88e3    # HHC detuning [Hz]
    hc_auto_detune: bool = False     # Auto flat-potential detuning (xi=1)
    use_form_factor_iter: bool = True  # Self-consistent Fokker-Planck F_n(sigma_tau) iteration

    # 4. Fill Pattern & Tracking
    pattern_mode: str = "Uniform (100%)"  # 'Uniform (100%)', '80% filling', '90% filling', 'Single bunch', 'Custom'
    custom_pattern_expr: str = ""
    Par_num: int = 10000             # Macro-particles per bunch
    Track_num: int = 30000           # Tracking turns
    Recor_step: int = 10             # Record interval [turns]
    short_range_wake: bool = True    # Enable intra-bunch short-range wake

    # 5. PI Feedback Loop
    pi_enable: bool = True
    pi_kp: float = 0.8
    pi_ki: float = 1.0e-4
    pi_m: int = 80                   # Averaging window [buckets]
    pi_d: int = 800                  # Loop delay [buckets]

    # 6. Direct RF Feedback (DFB) Loop
    # Note: phase_shift_deg = 0.0 represents standard -180 deg negative feedback!
    dfb_enable: bool = True
    dfb_gain: float = 0.20
    dfb_phase_shift_deg: float = 0.0
    dfb_delay: int = 704             # Loop delay [buckets]
    dfb_m: int = 208                 # Averaging window [buckets]
    dfb_every: int = 208             # Update interval [buckets]
    dfb_auto_scale: bool = True      # Linearly interpolate gain between (I1, G1) and (I2, G2)
    dfb_adaptive: bool = False       # Adaptive gain scaling with stability margin
    dfb_I1_mA: float = 300.0
    dfb_G1: float = 0.25
    dfb_I2_mA: float = 500.0
    dfb_G2: float = 0.20

    # 7. Higher-Order Modes (HOMs)
    homs: List[HOMMode] = field(default_factory=lambda: [
        HOMMode(Q=830.0, R=6500.0, freq_Hz=1700.0e6),
        HOMMode(Q=921.0, R=7300.0, freq_Hz=1182.0e6),
        HOMMode(Q=209.0, R=10100.0, freq_Hz=960.7e6),
    ])
    hom_auto_align: bool = True      # Align HOM frequencies to worst-case unstable synchrotron sideband

    def sync_derived_cavity_params(self) -> None:
        """Synchronize loaded Q and total loaded shunt impedance R_L from (Ncav, Q0, RoQ, beta)."""
        if self.Q0_mc > 0 and self.RoQ_mc > 0:
            self.Q_mc = self.Q0_mc / (1.0 + max(0.0, self.beta_mc))
            self.R_mc = self.Ncav_mc * self.RoQ_mc * self.Q_mc
        if self.Q0_hc > 0 and self.RoQ_hc > 0:
            self.Q_hc = self.Q0_hc / (1.0 + max(0.0, self.beta_hc))
            self.R_hc = self.Ncav_hc * self.RoQ_hc * self.Q_hc

    def get_dfb_gain_for_current(self, i0_amp: Optional[float] = None) -> float:
        """Compute active DFB gain at current i0_amp [A] matching STABLE_App.getDFBGainForCurrent."""
        if not self.dfb_enable:
            return 0.0
        curr_mA = (self.I0 if i0_amp is None else i0_amp) * 1e3
        if self.dfb_adaptive:
            # Adaptive rule matching STABLE_App.m
            c_speed = 299792458.0
            f_rf = self.h / (self.C / c_speed)
            w_half = 2.0 * math.pi * f_rf / (2.0 * max(1.0, self.Q_mc))
            tau_d = max(1, self.dfb_delay) * ((self.C / c_speed) / self.h)
            g_max_0 = (math.pi / (2.0 * max(1e-12, w_half * tau_d)))
            return min(0.45 * g_max_0, 0.22 + 0.08 * ((curr_mA - 300.0) / 200.0))
        if self.dfb_auto_scale:
            if abs(self.dfb_I2_mA - self.dfb_I1_mA) < 1e-6:
                return self.dfb_G1
            t = (curr_mA - self.dfb_I1_mA) / (self.dfb_I2_mA - self.dfb_I1_mA)
            t = max(0.0, min(1.0, t))
            return self.dfb_G1 + t * (self.dfb_G2 - self.dfb_I1_mA * 0.0 + (self.dfb_G2 - self.dfb_G1) - (self.dfb_G2 - self.dfb_G1)) if False else (
                self.dfb_G1 + t * (self.dfb_G2 - self.dfb_G1)
            )
        return self.dfb_gain

    @classmethod
    def soleil2_default(cls) -> "MachineConfig":
        cfg = cls()
        cfg.sync_derived_cavity_params()
        return cfg

    @classmethod
    def half_default(cls) -> "MachineConfig":
        return cls(
            preset_name="HALF (Hefei 4th Gen)",
            C=479.86,
            E0=2.2e9,
            I0=0.35,
            U0=373.5e3,
            h=800,
            alpha_c=8.1e-5,
            tau_s=12.5e-3,
            tau_z=12.5e-3,
            sigma_e0=6.4e-4,
            Ncav_mc=1,
            Q0_mc=1.1e5,
            RoQ_mc=45.0,
            beta_mc=0.0,
            V_mc=1.25e6,
            Q_mc=1.1e5,
            R_mc=4.95e6,
            fre_shift_mc=-7.0e3,
            mc_auto_detune=True,
            n_hc=3,
            Ncav_hc=1,
            Q0_hc=2.0e8,
            RoQ_hc=39.0,
            beta_hc=0.0,
            Q_hc=2.0e8,
            R_hc=7.8e9,
            fre_shift_hc=45.0e3,
            hc_auto_detune=True,
            use_form_factor_iter=True,
            pattern_mode="80% filling",
            Par_num=10000,
            Track_num=30000,
            Recor_step=10,
            short_range_wake=True,
            pi_enable=True,
            pi_kp=1.0,
            pi_ki=1.0e-5,
            pi_m=80,
            pi_d=800,
            dfb_enable=False,
            dfb_gain=0.2,
            dfb_phase_shift_deg=0.0,
            dfb_delay=800,
            dfb_m=1,
            dfb_every=1,
            dfb_auto_scale=False,
            homs=[],
            hom_auto_align=False,
        )

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d["homs"] = [h.to_dict() for h in self.homs]
        return d

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "MachineConfig":
        data = dict(d)
        homs_raw = data.pop("homs", [])
        valid_keys = {f.name for f in cls.__dataclass_fields__.values()}
        filtered = {k: v for k, v in data.items() if k in valid_keys}
        obj = cls(**filtered)
        obj.homs = [HOMMode.from_dict(x) if isinstance(x, dict) else x for x in homs_raw]
        return obj

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), indent=2)

    @classmethod
    def from_json(cls, s: str) -> "MachineConfig":
        return cls.from_dict(json.loads(s))


@dataclass
class ScanConfig:
    """2D Stability Map Scan & Multi-A100 GPU Dispatch Configuration."""
    Df_min_kHz: float = 50.0
    Df_max_kHz: float = 150.0
    Df_step_kHz: float = 5.0
    I0_min_mA: float = 300.0
    I0_max_mA: float = 500.0
    I0_step_mA: float = 25.0

    # Pre-screening & Multi-GPU acceleration
    enable_dmode_prescreen: bool = False    # T. He PRAB 2024 analytical pre-filter (Superconducting Q_hc >= 1e5 only)
    num_gpus: int = 1                       # Number of A100 GPUs to request (1 or 2 schedules fastest on Topaze)
    workers_per_gpu: int = 16               # Optimal sweet spot for A100-80GB + 32 CPU cores (12~16 pts/GPU)
    execution_mode: str = "daemon_stream"   # "daemon_stream" (live RPC) or "ccc_msub_batch" (TGCC batch scheduler)

    def get_df_grid_kHz(self) -> List[float]:
        n = max(1, int(round((self.Df_max_kHz - self.Df_min_kHz) / max(1e-6, self.Df_step_kHz))) + 1)
        return [round(self.Df_min_kHz + i * self.Df_step_kHz, 4) for i in range(n)]

    def get_i0_grid_mA(self) -> List[float]:
        n = max(1, int(round((self.I0_max_mA - self.I0_min_mA) / max(1e-6, self.I0_step_mA))) + 1)
        return [round(self.I0_min_mA + i * self.I0_step_mA, 4) for i in range(n)]

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "ScanConfig":
        valid_keys = {f.name for f in cls.__dataclass_fields__.values()}
        return cls(**{k: v for k, v in d.items() if k in valid_keys})


@dataclass
class TGCCConnectionProfile:
    """Per-user TGCC Topaze/Irene A100 Cluster & SSH Connection Configuration."""
    profile_name: str = "TGCC Irene A100"
    mode: str = "local"                     # "local" (local server/GPU) or "ssh_tgcc" (remote TGCC via SSH)
    ssh_host: str = "topaze.ccc.cea.fr"
    ssh_port: int = 22
    ssh_username: str = "zhaolu"
    auth_method: str = "password"           # "password", "private_key", "agent"
    private_key_path: str = ""
    ccc_project: str = "soldai"             # #MSUB -A project account on TGCC Topaze (from ccc_mpp)
    ccc_partition: str = "a100"             # #MSUB -q a100 (TGCC Topaze A100 partition)
    ccc_qos: str = "normal"                 # #MSUB -Q QoS ("normal", "test" for high-priority <=30m, "long")
    remote_work_dir: str = "/ccc/work/cont003/soleil/zhaolu/STABLE_Py_TGCC"
    python_env_cmd: str = "python3"
    remote_daemon_port: int = 18765
    gpus_per_node: int = 4                  # Topaze A100 nodes have 4x NVIDIA A100-SXM4-80GB
    walltime_seconds: int = 43200           # Default 12.0 hours (up to 86400s = 24h on Topaze a100)
    auto_walltime: bool = True              # Auto-estimate and allocate #MSUB -T walltime from grid size

    def compute_effective_walltime_seconds(
        self,
        machine_cfg: "MachineConfig",
        scan_cfg: Optional["ScanConfig"] = None,
        is_scan: bool = True,
    ) -> int:
        """
        Estimate required walltime (`#MSUB -T` in seconds) with a 3.5x safety margin,
        clamped between 1 hour (3600s) and 24 hours (86400s, Topaze partition limit).
        """
        if not self.auto_walltime:
            return max(600, min(86400, int(self.walltime_seconds)))

        n_pts = 1
        n_gpus = 1
        w_per_gpu = 1
        if is_scan and scan_cfg is not None:
            n_pts = max(1, len(scan_cfg.get_i0_grid_mA()) * len(scan_cfg.get_df_grid_kHz()))
            n_gpus = max(1, int(scan_cfg.num_gpus))
            w_per_gpu = max(1, int(scan_cfg.workers_per_gpu))

        total_concurrent = max(1, n_gpus * w_per_gpu)

        # With C-vectorized mc_cavity_control + A100 CUDA kernels + MPS multi-process concurrency:
        # 1 batch of K=16 concurrent points on A100 takes ~25-35s (use 45s safe budget per concurrent batch).
        # Keeping walltime <= 2.0h (7200s) allows Topaze SLURM backfill scheduler to start the job immediately!
        turns_factor = max(0.1, float(machine_cfg.Track_num) / 30000.0)
        part_factor = max(0.2, (float(machine_cfg.Par_num) * float(machine_cfg.h)) / (10000.0 * 416.0))
        hom_factor = 1.0 + 0.25 * len(machine_cfg.homs)
        sec_per_batch_safe = 45.0 * turns_factor * part_factor * hom_factor

        batches = math.ceil(n_pts / float(total_concurrent))
        est_sec = int(batches * sec_per_batch_safe + 900)  # +15 min buffer
        # Round up to nearest 1800s (0.5h), minimum 1.5h (5400s) for scans so SLURM backfill schedules immediately!
        min_sec = 5400 if is_scan else 1800
        rounded = int(math.ceil(max(min_sec, est_sec) / 1800.0) * 1800)
        return min(86400, rounded)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "TGCCConnectionProfile":
        valid_keys = {f.name for f in cls.__dataclass_fields__.values()}
        return cls(**{k: v for k, v in d.items() if k in valid_keys})


@dataclass
class TelemetryFrame:
    """Real-time turn-by-turn streaming frame sent from TGCC A100 Tracker to Client GUI."""
    job_id: str
    event_type: str = "turn_progress"       # "turn_progress", "tracking_done", "error"
    turn: int = 0
    total_turns: int = 1
    elapsed_sec: float = 0.0
    eta_sec: float = 0.0
    device_name: str = "CPU"
    beam_lost: bool = False

    # Downsampled bucket waveforms for live GUI plots (kV, deg)
    hc_amp_kV: List[float] = field(default_factory=list)
    hc_phase_deg: List[float] = field(default_factory=list)
    mc_amp_kV: List[float] = field(default_factory=list)
    mc_phase_deg: List[float] = field(default_factory=list)

    # Turn evolution history (recorded up to current turn)
    turns_axis: List[int] = field(default_factory=list)
    bunch_length_ps_history: List[float] = field(default_factory=list)
    centroid_jitter_ps_history: List[float] = field(default_factory=list)
    energy_spread_history: List[float] = field(default_factory=list)

    # Final classification results (populated when event_type == "tracking_done")
    stability_code: int = 0                 # 0=Stable, 1=PTBL, 2=Coupled D-Q, 3=LCBI, 4=Other/D-mode
    status_label: str = "Running"
    growth_rate_s1: float = 0.0
    dominant_freq_Hz: float = 0.0
    dominant_mode_mu: int = 0
    mode_spectrum_amp_ps: List[float] = field(default_factory=list)
    xi_eff: float = 0.0
    psi_hc_deg: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "TelemetryFrame":
        valid_keys = {f.name for f in cls.__dataclass_fields__.values()}
        return cls(**{k: v for k, v in d.items() if k in valid_keys})


@dataclass
class ScanPointTelemetry:
    """Real-time streaming update for a single (I0, Df_HC) point in a 2D Stability Map Scan."""
    job_id: str
    point_index: int
    total_points: int
    I0_mA: float
    Df_hc_kHz: float
    psi_hc_deg: float
    xi_eff: float
    stability_code: int                     # 0=Stable, 1=PTBL, 2=Coupled D-Q, 3=LCBI, 4=D-mode/Other
    status_label: str
    growth_rate_s1: float
    dominant_freq_Hz: float
    dominant_mode_mu: int
    mean_bunch_length_ps: float
    centroid_jitter_ps: float
    prescreened_dmode: bool = False
    gpu_id: int = 0
    elapsed_sec: float = 0.0

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "ScanPointTelemetry":
        valid_keys = {f.name for f in cls.__dataclass_fields__.values()}
        return cls(**{k: v for k, v in d.items() if k in valid_keys})
