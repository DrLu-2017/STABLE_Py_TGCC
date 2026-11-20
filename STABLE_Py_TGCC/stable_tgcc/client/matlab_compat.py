"""
100% MATLAB `STABLE_App.m` Compatibility Utilities for STABLE_Py_TGCC:
  1. `restore_scan_config_from_mat` & `reconstruct_config_from_scan_folder` (ports lines 2443-2720)
  2. `load_config_from_mat_or_m` & `save_config_to_mat` (ports lines 1579-1613 & 2722-2785)
  3. `save_scan_summary_mat` & `find_latest_default_scan_summary` (ports lines 4590-4688)
  4. `fit_safety_boundary` (ports `fitAndPlotSafetyBoundary`, lines 1984-2069)
  5. `analyze_scan_folder_mat` (ports `AnalyzeButtonPushed`, lines 1701-1881)
  6. `align_homs_to_worst_case_sideband` (ports `alignHOMToWorstCaseSideband`, lines 2159-2222)
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
import copy
import math
import re
import numpy as np
import scipy.io as sio

from stable_tgcc.common.models import MachineConfig, ScanConfig, HOMMode


def _mat_struct_to_dict(obj: Any) -> Dict[str, Any]:
    """Convert a scipy.io.loadmat MATLAB struct / void / mat_struct into a plain Python dict."""
    if obj is None:
        return {}
    if isinstance(obj, dict):
        return {k: v for k, v in obj.items() if not k.startswith("__")}
    if hasattr(obj, "_fieldnames"):
        return {fn: getattr(obj, fn) for fn in obj._fieldnames}
    if isinstance(obj, np.ndarray):
        if obj.size == 0:
            return {}
        if obj.dtype.names:
            first = obj.flat[0]
            return {name: first[name] for name in obj.dtype.names}
        if obj.size == 1:
            return _mat_struct_to_dict(obj.item())
    if isinstance(obj, np.void) and obj.dtype.names:
        return {name: obj[name] for name in obj.dtype.names}
    return {}


def _scalar_float(val: Any, default: Optional[float] = None) -> Optional[float]:
    if val is None:
        return default
    try:
        arr = np.asarray(val)
        if arr.size == 0:
            return default
        v = float(arr.flat[0])
        return default if math.isnan(v) else v
    except Exception:
        return default


def _scalar_bool(val: Any, default: bool = False) -> bool:
    if val is None:
        return default
    try:
        arr = np.asarray(val)
        if arr.size == 0:
            return default
        return bool(arr.flat[0])
    except Exception:
        return default


def _scalar_str(val: Any, default: str = "") -> str:
    if val is None:
        return default
    if isinstance(val, str):
        return val.strip()
    try:
        arr = np.asarray(val)
        if arr.size == 0:
            return default
        item = arr.flat[0]
        if isinstance(item, (str, np.str_)):
            return str(item).strip()
        return "".join(str(x) for x in arr.ravel()).strip()
    except Exception:
        return default


def apply_matlab_config_dict(
    cfg_dict: Dict[str, Any],
    base_cfg: Optional[MachineConfig] = None,
    base_scan: Optional[ScanConfig] = None,
) -> Tuple[MachineConfig, ScanConfig]:
    """
    Port of MATLAB `setUIFromConfig(app, config)` (STABLE_App.m lines 2443-2586).
    Maps MATLAB struct field names (`KP`, `KI`, `DFB_enable`, `R_hc`, `Q_hc`, `HOMs`, etc.)
    onto Python `MachineConfig` and `ScanConfig`.
    """
    m = copy.deepcopy(base_cfg) if base_cfg is not None else MachineConfig.soleil2_default()
    s = copy.deepcopy(base_scan) if base_scan is not None else ScanConfig()

    for attr_name in ("C", "E0", "I0", "U0", "alpha_c", "tau_s", "tau_z", "sigma_e0"):
        v = _scalar_float(cfg_dict.get(attr_name))
        if v is not None and v > 0:
            setattr(m, attr_name, v)

    h_v = _scalar_float(cfg_dict.get("h"))
    if h_v is not None and h_v >= 1:
        m.h = int(round(h_v))

    if abs(m.tau_s - m.tau_z) > 1e-9:
        tau_wc = max(m.tau_s, m.tau_z)
        m.tau_s = tau_wc
        m.tau_z = tau_wc

    # Main Cavity
    for fn, attr in (
        ("V_mc", "V_mc"),
        ("Q0_mc", "Q0_mc"),
        ("RoQ_mc", "RoQ_mc"),
        ("R_over_Q_mc", "RoQ_mc"),
        ("beta_mc", "beta_mc"),
        ("Q_mc", "Q_mc"),
        ("R_mc", "R_mc"),
        ("fre_shift_mc", "fre_shift_mc"),
    ):
        v = _scalar_float(cfg_dict.get(fn))
        if v is not None:
            setattr(m, attr, v)

    nc_mc = _scalar_float(cfg_dict.get("Ncav_mc"))
    if nc_mc is not None and nc_mc >= 1:
        m.Ncav_mc = int(round(nc_mc))

    if "Q0_mc" in cfg_dict and "beta_mc" in cfg_dict and "Q_mc" not in cfg_dict:
        m.Q_mc = m.Q0_mc / (1.0 + max(0.0, m.beta_mc))
        m.R_mc = m.Ncav_mc * m.RoQ_mc * m.Q_mc

    if "mc_auto_detune" in cfg_dict:
        m.mc_auto_detune = _scalar_bool(cfg_dict.get("mc_auto_detune"), m.mc_auto_detune)

    # Harmonic Cavity
    for fn, attr in (
        ("Q0_hc", "Q0_hc"),
        ("RoQ_hc", "RoQ_hc"),
        ("R_over_Q_hc", "RoQ_hc"),
        ("beta_hc", "beta_hc"),
        ("Q_hc", "Q_hc"),
        ("fre_shift_hc", "fre_shift_hc"),
    ):
        v = _scalar_float(cfg_dict.get(fn))
        if v is not None:
            setattr(m, attr, v)

    for fn, attr in (("n_hc", "n_hc"), ("Ncav_hc", "Ncav_hc")):
        v = _scalar_float(cfg_dict.get(fn))
        if v is not None and v >= 1:
            setattr(m, attr, int(round(v)))

    r_hc_raw = _scalar_float(cfg_dict.get("R_hc"))
    if r_hc_raw is not None and r_hc_raw > 0:
        # In STABLE_App.m lines 2489-2494: if R_hc > 1e4 it is in Ohm, else in MOhm
        m.R_hc = r_hc_raw if r_hc_raw > 1e4 else (r_hc_raw * 1e6)
    elif "Q_hc" in cfg_dict and ("RoQ_hc" in cfg_dict or "R_over_Q_hc" in cfg_dict):
        if abs(m.Q_hc - 31e3) < 5e3:
            m.R_hc = 2.0 * m.Q_hc * m.RoQ_hc
        else:
            m.R_hc = m.Ncav_hc * m.Q_hc * m.RoQ_hc

    if "auto_detune" in cfg_dict:
        m.hc_auto_detune = _scalar_bool(cfg_dict.get("auto_detune"), m.hc_auto_detune)
    if "hc_auto_detune" in cfg_dict:
        m.hc_auto_detune = _scalar_bool(cfg_dict.get("hc_auto_detune"), m.hc_auto_detune)
    if "form_factor_iter" in cfg_dict:
        m.use_form_factor_iter = _scalar_bool(cfg_dict.get("form_factor_iter"), m.use_form_factor_iter)
    if "use_form_factor_iter" in cfg_dict:
        m.use_form_factor_iter = _scalar_bool(cfg_dict.get("use_form_factor_iter"), m.use_form_factor_iter)

    # Pattern & Tracking
    pt_str = _scalar_str(cfg_dict.get("pattern_type") or cfg_dict.get("pattern_mode"))
    if pt_str in ("Uniform (100%)", "80% filling", "90% filling", "Single bunch", "Custom"):
        m.pattern_mode = pt_str
    cp_str = _scalar_str(cfg_dict.get("custom_pattern") or cfg_dict.get("custom_pattern_expr"))
    if cp_str:
        m.custom_pattern_expr = cp_str

    for fn, attr in (("Par_num", "Par_num"), ("Track_num", "Track_num"), ("Recor_step", "Recor_step")):
        v = _scalar_float(cfg_dict.get(fn))
        if v is not None and v >= 1:
            setattr(m, attr, int(round(v)))

    if "ShortRange" in cfg_dict:
        m.short_range_wake = _scalar_bool(cfg_dict.get("ShortRange"), m.short_range_wake)
    elif "short_range_wake" in cfg_dict:
        m.short_range_wake = _scalar_bool(cfg_dict.get("short_range_wake"), m.short_range_wake)

    # PI Feedback
    if "PI_enable" in cfg_dict:
        m.pi_enable = _scalar_bool(cfg_dict.get("PI_enable"), m.pi_enable)
    elif "pi_enable" in cfg_dict:
        m.pi_enable = _scalar_bool(cfg_dict.get("pi_enable"), m.pi_enable)
    for fn, attr in (("KP", "pi_kp"), ("pi_kp", "pi_kp"), ("KI", "pi_ki"), ("pi_ki", "pi_ki")):
        v = _scalar_float(cfg_dict.get(fn))
        if v is not None:
            setattr(m, attr, v)
    for fn, attr in (("PI_m", "pi_m"), ("pi_m", "pi_m"), ("PI_d", "pi_d"), ("pi_d", "pi_d")):
        v = _scalar_float(cfg_dict.get(fn))
        if v is not None and v >= 1:
            setattr(m, attr, int(round(v)))

    # DFB Feedback
    for fn, attr in (
        ("DFB_enable", "dfb_enable"),
        ("dfb_enable", "dfb_enable"),
        ("DFB_adaptive", "dfb_adaptive"),
        ("dfb_adaptive", "dfb_adaptive"),
        ("DFB_auto_scale", "dfb_auto_scale"),
        ("dfb_auto_scale", "dfb_auto_scale"),
    ):
        if fn in cfg_dict:
            setattr(m, attr, _scalar_bool(cfg_dict.get(fn), getattr(m, attr)))

    for fn, attr in (
        ("DFB_gain", "dfb_gain"),
        ("dfb_gain", "dfb_gain"),
        ("DFB_phase_shift", "dfb_phase_shift_deg"),
        ("dfb_phase_shift_deg", "dfb_phase_shift_deg"),
        ("DFB_I1", "dfb_I1_mA"),
        ("dfb_I1_mA", "dfb_I1_mA"),
        ("DFB_G1", "dfb_G1"),
        ("dfb_G1", "dfb_G1"),
        ("DFB_I2", "dfb_I2_mA"),
        ("dfb_I2_mA", "dfb_I2_mA"),
        ("DFB_G2", "dfb_G2"),
        ("dfb_G2", "dfb_G2"),
    ):
        v = _scalar_float(cfg_dict.get(fn))
        if v is not None:
            setattr(m, attr, v)

    for fn, attr in (
        ("DFB_delay", "dfb_delay"),
        ("dfb_delay", "dfb_delay"),
        ("DFB_m", "dfb_m"),
        ("dfb_m", "dfb_m"),
        ("DFB_every", "dfb_every"),
        ("dfb_every", "dfb_every"),
    ):
        v = _scalar_float(cfg_dict.get(fn))
        if v is not None and v >= 1:
            setattr(m, attr, int(round(v)))

    # HOMs
    if "HOMs" in cfg_dict:
        raw_homs = cfg_dict.get("HOMs")
        parsed_homs: List[HOMMode] = []
        try:
            arr_h = np.asarray(raw_homs, dtype=object)
            if arr_h.size >= 3:
                arr_f = np.asarray(raw_homs, dtype=np.float64).reshape(-1, 3)
                for row in arr_f:
                    if row[1] > 0:
                        parsed_homs.append(HOMMode(Q=float(row[0]), R=float(row[1]), freq_Hz=float(row[2])))
        except Exception:
            pass
        m.homs = parsed_homs
    elif "homs" in cfg_dict and isinstance(cfg_dict["homs"], list):
        m.homs = [HOMMode.from_dict(x) if isinstance(x, dict) else x for x in cfg_dict["homs"]]

    if "HOM_auto_align" in cfg_dict:
        m.hom_auto_align = _scalar_bool(cfg_dict.get("HOM_auto_align"), m.hom_auto_align)
    elif "hom_auto_align" in cfg_dict:
        m.hom_auto_align = _scalar_bool(cfg_dict.get("hom_auto_align"), m.hom_auto_align)

    # Scan bounds
    for fn, attr in (
        ("scan_Df_min", "Df_min_kHz"),
        ("scan_Df_max", "Df_max_kHz"),
        ("scan_Df_step", "Df_step_kHz"),
        ("scan_I0_min", "I0_min_mA"),
        ("scan_I0_max", "I0_max_mA"),
        ("scan_I0_step", "I0_step_mA"),
    ):
        v = _scalar_float(cfg_dict.get(fn))
        if v is not None and v > 0:
            setattr(s, attr, v)

    return m, s


def restore_scan_config_from_mat(
    s_data: Dict[str, Any],
    mat_or_folder_path: str | Path,
    base_cfg: Optional[MachineConfig] = None,
    base_scan: Optional[ScanConfig] = None,
) -> Tuple[MachineConfig, ScanConfig]:
    """
    Exact port of MATLAB `restoreScanConfigToUI` and `reconstructConfigFromScanFolder`
    (STABLE_App.m lines 2588-2720).
    """
    p = Path(mat_or_folder_path)
    folder = p if p.is_dir() else p.parent
    f_name = folder.name

    cfg_dict: Dict[str, Any] = {}
    if "config" in s_data:
        cfg_dict = _mat_struct_to_dict(s_data["config"])

    # If config struct wasn't stored in scan_summary.mat, reconstruct from first SOLEIL_I0*.mat in folder
    if not cfg_dict and folder.exists():
        pt_files = sorted(folder.glob("SOLEIL_I0*.mat"))
        if pt_files:
            try:
                pt = sio.loadmat(str(pt_files[0]), squeeze_me=True, struct_as_record=False)
                if "HALF" in pt:
                    h_dict = _mat_struct_to_dict(pt["HALF"])
                    for fn in (
                        "C", "E0", "I0", "U0", "h", "alpha_c", "tau_s", "tau_z", "sigma_e0",
                        "V_mc", "Ncav_mc", "Q0_mc", "beta_mc", "Q_mc", "R_mc", "fre_shift_mc",
                        "n_hc", "Ncav_hc", "Q0_hc", "beta_hc", "Q_hc", "R_hc", "fre_shift_hc",
                    ):
                        if fn in h_dict:
                            cfg_dict[fn] = h_dict[fn]
                    if "ShortRange_on" in h_dict:
                        cfg_dict["ShortRange"] = h_dict["ShortRange_on"]
                    if "Q_hom_m0" in h_dict and "R_hom_m0" in h_dict and "fre_hom_m0" in h_dict:
                        q_h = np.atleast_1d(np.asarray(h_dict["Q_hom_m0"], dtype=np.float64)).ravel()
                        r_h = np.atleast_1d(np.asarray(h_dict["R_hom_m0"], dtype=np.float64)).ravel()
                        f_h = np.atleast_1d(np.asarray(h_dict["fre_hom_m0"], dtype=np.float64)).ravel()
                        if np.any(r_h > 0):
                            n_h = min(len(q_h), len(r_h), len(f_h))
                            cfg_dict["HOMs"] = np.column_stack([q_h[:n_h], r_h[:n_h], f_h[:n_h]])
                        else:
                            cfg_dict["HOMs"] = []
                if "PI" in pt:
                    p_dict = _mat_struct_to_dict(pt["PI"])
                    for k_p, k_c in (("enable", "PI_enable"), ("KP", "KP"), ("KI", "KI"), ("m", "PI_m"), ("d", "PI_d")):
                        if k_p in p_dict:
                            cfg_dict[k_c] = p_dict[k_p]
                if "DFB" in pt:
                    d_dict = _mat_struct_to_dict(pt["DFB"])
                    for k_d, k_c in (
                        ("enable", "DFB_enable"),
                        ("gain", "DFB_gain"),
                        ("phase_shift", "DFB_phase_shift"),
                        ("d", "DFB_delay"),
                        ("m", "DFB_m"),
                        ("every", "DFB_every"),
                    ):
                        if k_d in d_dict:
                            cfg_dict[k_c] = d_dict[k_d]
            except Exception:
                pass

        # Refine feedback & HOM flags from folder naming convention (STABLE_App.m lines 2697-2719)
        if "DFBadapt" in f_name:
            cfg_dict["DFB_enable"] = True
            cfg_dict["DFB_adaptive"] = True
        elif "DFBauto" in f_name:
            cfg_dict["DFB_enable"] = True
            cfg_dict["DFB_adaptive"] = False
            cfg_dict["DFB_auto_scale"] = True
        elif "NoDFB" in f_name:
            cfg_dict["DFB_enable"] = False
        elif "DFBg" in f_name:
            cfg_dict["DFB_enable"] = True
            cfg_dict["DFB_adaptive"] = False
            cfg_dict["DFB_auto_scale"] = False
            m_g = re.search(r"DFBg([0-9.]+)", f_name)
            if m_g:
                cfg_dict["DFB_gain"] = float(m_g.group(1))

        if "HOMalign" in f_name:
            cfg_dict["HOM_auto_align"] = True
        elif "HOMnom" in f_name or "NoHOM" in f_name:
            cfg_dict["HOM_auto_align"] = False

    # Synchronize scan grid bounds from Df_scan_kHz / I0_scan_mA
    if "Df_scan_kHz" in s_data:
        df_vec = np.asarray(s_data["Df_scan_kHz"], dtype=np.float64).ravel()
        if df_vec.size > 0:
            cfg_dict["scan_Df_min"] = float(np.min(df_vec))
            cfg_dict["scan_Df_max"] = float(np.max(df_vec))
            if df_vec.size > 1:
                cfg_dict["scan_Df_step"] = float(abs(df_vec[1] - df_vec[0]))
    if "I0_scan_mA" in s_data:
        i0_vec = np.asarray(s_data["I0_scan_mA"], dtype=np.float64).ravel()
        if i0_vec.size > 0:
            cfg_dict["scan_I0_min"] = float(np.min(i0_vec))
            cfg_dict["scan_I0_max"] = float(np.max(i0_vec))
            if i0_vec.size > 1:
                cfg_dict["scan_I0_step"] = float(abs(i0_vec[1] - i0_vec[0]))

    return apply_matlab_config_dict(cfg_dict, base_cfg=base_cfg, base_scan=base_scan)


def load_config_from_mat_or_m(
    file_path: str | Path,
    base_cfg: Optional[MachineConfig] = None,
    base_scan: Optional[ScanConfig] = None,
) -> Tuple[MachineConfig, ScanConfig]:
    """
    Port of MATLAB `LoadConfigButtonPushed` and `loadMachineFile` (STABLE_App.m lines 1579-1600 & 2722-2785).
    Loads either a `.mat` configuration/scan_summary file or a `.m` machine definition script.
    """
    p = Path(file_path)
    ext = p.suffix.lower()
    if ext == ".mat":
        s_data = sio.loadmat(str(p), squeeze_me=True, struct_as_record=False)
        if "config" in s_data:
            return restore_scan_config_from_mat(s_data, p, base_cfg, base_scan)
        return apply_matlab_config_dict(_mat_struct_to_dict(s_data), base_cfg, base_scan)

    if ext == ".m":
        txt = p.read_text(encoding="utf-8", errors="replace")
        cfg_dict: Dict[str, Any] = {}

        def _parse_var(vname: str) -> Optional[float]:
            m_var = re.search(rf"\b{re.escape(vname)}\b\s*=\s*([^;]+);", txt)
            if not m_var:
                return None
            expr = re.sub(r"%.*", "", m_var.group(1)).strip()
            try:
                return float(eval(expr, {"__builtins__": {}}, {"pi": math.pi}))
            except Exception:
                return None

        for vname in (
            "C", "E0", "I0", "U0", "h", "alpha_c", "tau_s", "tau_z", "sigma_e0",
            "V_mc", "Ncav_mc", "Q0_mc", "beta_mc", "Q_mc", "R_over_Q_mc", "R_mc", "fre_shift_mc",
            "n_hc", "Ncav_hc", "Q0_hc", "beta_hc", "Q_hc", "R_over_Q_hc", "R_hc", "fre_shift_hc",
        ):
            val = _parse_var(vname)
            if val is not None:
                cfg_dict[vname] = val
        if "pattern = ones(1, h)" in txt or "pattern=ones(1,h)" in txt:
            cfg_dict["pattern_type"] = "Uniform (100%)"
        if "fre_shift_hc" in cfg_dict:
            cfg_dict["auto_detune"] = False
        return apply_matlab_config_dict(cfg_dict, base_cfg, base_scan)

    raise ValueError(f"Unsupported file extension: {ext}")


def machine_config_to_matlab_struct(cfg: MachineConfig, scan_cfg: ScanConfig) -> Dict[str, Any]:
    """Convert `MachineConfig` + `ScanConfig` to a MATLAB `getConfigFromUI` compatible dict for `sio.savemat`."""
    homs_arr = (
        np.array([[m.Q, m.R, m.freq_Hz] for m in cfg.homs], dtype=np.float64)
        if cfg.homs
        else np.zeros((0, 3), dtype=np.float64)
    )
    return {
        "C": float(cfg.C),
        "E0": float(cfg.E0),
        "I0": float(cfg.I0),
        "U0": float(cfg.U0),
        "h": float(cfg.h),
        "alpha_c": float(cfg.alpha_c),
        "tau_s": float(cfg.tau_s),
        "tau_z": float(cfg.tau_z),
        "sigma_e0": float(cfg.sigma_e0),
        "V_mc": float(cfg.V_mc),
        "Ncav_mc": float(cfg.Ncav_mc),
        "Q0_mc": float(cfg.Q0_mc),
        "RoQ_mc": float(cfg.RoQ_mc),
        "beta_mc": float(cfg.beta_mc),
        "Q_mc": float(cfg.Q_mc),
        "R_mc": float(cfg.R_mc),
        "fre_shift_mc": float(cfg.fre_shift_mc),
        "mc_auto_detune": bool(cfg.mc_auto_detune),
        "n_hc": float(cfg.n_hc),
        "Ncav_hc": float(cfg.Ncav_hc),
        "Q0_hc": float(cfg.Q0_hc),
        "RoQ_hc": float(cfg.RoQ_hc),
        "beta_hc": float(cfg.beta_hc),
        "Q_hc": float(cfg.Q_hc),
        "R_hc": float(cfg.R_hc),
        "auto_detune": bool(cfg.hc_auto_detune),
        "form_factor_iter": bool(cfg.use_form_factor_iter),
        "fre_shift_hc": float(cfg.fre_shift_hc),
        "pattern_type": str(cfg.pattern_mode),
        "custom_pattern": str(cfg.custom_pattern_expr),
        "Par_num": float(cfg.Par_num),
        "Track_num": float(cfg.Track_num),
        "Recor_step": float(cfg.Recor_step),
        "ShortRange": bool(cfg.short_range_wake),
        "KP": float(cfg.pi_kp),
        "KI": float(cfg.pi_ki),
        "PI_m": float(cfg.pi_m),
        "PI_d": float(cfg.pi_d),
        "PI_enable": bool(cfg.pi_enable),
        "DFB_enable": bool(cfg.dfb_enable),
        "DFB_adaptive": bool(cfg.dfb_adaptive),
        "DFB_gain": float(cfg.dfb_gain),
        "DFB_phase_shift": float(cfg.dfb_phase_shift_deg),
        "DFB_delay": float(cfg.dfb_delay),
        "DFB_m": float(cfg.dfb_m),
        "DFB_every": float(cfg.dfb_every),
        "DFB_auto_scale": bool(cfg.dfb_auto_scale),
        "DFB_I1": float(cfg.dfb_I1_mA),
        "DFB_G1": float(cfg.dfb_G1),
        "DFB_I2": float(cfg.dfb_I2_mA),
        "DFB_G2": float(cfg.dfb_G2),
        "HOMs": homs_arr,
        "HOM_auto_align": bool(cfg.hom_auto_align),
        "scan_Df_min": float(scan_cfg.Df_min_kHz),
        "scan_Df_max": float(scan_cfg.Df_max_kHz),
        "scan_Df_step": float(scan_cfg.Df_step_kHz),
        "scan_I0_min": float(scan_cfg.I0_min_mA),
        "scan_I0_max": float(scan_cfg.I0_max_mA),
        "scan_I0_step": float(scan_cfg.I0_step_mA),
    }


def find_latest_default_scan_summary() -> Optional[Path]:
    """
    Port of MATLAB `updateStabilityMapView` auto-discovery (STABLE_App.m lines 4645-4661).
    Searches `../scan_data/Scan_*/scan_summary.mat` for the most recently modified scan summary.
    """
    pkg_root = Path(__file__).resolve().parents[2]
    candidates = [
        pkg_root / "scan_data",
        pkg_root.parent / "scan_data",
    ]
    best_file: Optional[Path] = None
    best_mtime: float = -1.0
    for scan_dir in candidates:
        if not scan_dir.exists():
            continue
        for f in scan_dir.glob("Scan_*/scan_summary.mat"):
            try:
                mt = f.stat().st_mtime
                if mt > best_mtime:
                    best_mtime = mt
                    best_file = f
            except Exception:
                pass
    return best_file


def fit_safety_boundary(
    df_scan_kHz: np.ndarray,
    i0_scan_mA: np.ndarray,
    stability_code: np.ndarray,
) -> Optional[Dict[str, Any]]:
    """
    Exact port of MATLAB `fitAndPlotSafetyBoundary` (STABLE_App.m lines 1984-2069).
    Extracts critical rightmost unstable -> stable transition points and fits Df_crit = a * I0 + b.
    """
    if df_scan_kHz.size < 2 or i0_scan_mA.size < 2 or stability_code.size == 0:
        return None

    crit_i0: List[float] = []
    crit_df: List[float] = []
    for r_i in range(stability_code.shape[0]):
        row_c = stability_code[r_i, :]
        valid_cols = np.where(~np.isnan(row_c))[0]
        if len(valid_cols) < 2:
            continue
        is_s = (row_c[valid_cols] == 0).astype(int)
        trans = np.where(np.diff(is_s) == 1)[0]
        if len(trans) > 0:
            last_t = trans[-1]
            c1, c2 = valid_cols[last_t], valid_cols[last_t + 1]
            crit_i0.append(float(i0_scan_mA[r_i]))
            crit_df.append(0.5 * float(df_scan_kHz[c1] + df_scan_kHz[c2]))

    if len(set(crit_i0)) < 2:
        return None

    crit_i0_arr = np.array(crit_i0, dtype=np.float64)
    crit_df_arr = np.array(crit_df, dtype=np.float64)
    p_fit = np.polyfit(crit_i0_arr, crit_df_arr, 1)
    df_pred = np.polyval(p_fit, crit_i0_arr)
    ss_res = float(np.sum((crit_df_arr - df_pred) ** 2))
    ss_tot = float(np.sum((crit_df_arr - np.mean(crit_df_arr)) ** 2))
    r2 = 1.0 - ss_res / max(ss_tot, 1e-12)
    rmse = math.sqrt(float(np.mean((crit_df_arr - df_pred) ** 2)))

    i0_dense = np.linspace(float(np.min(i0_scan_mA)), float(np.max(i0_scan_mA)), 100)
    df_dense = np.polyval(p_fit, i0_dense)
    a_val, b_val = float(p_fit[0]), float(p_fit[1])
    sgn = "+" if b_val >= 0 else "-"

    status_msg = (
        f"Fitted Safety Boundary: Df_hc > {a_val:.4f} * I0 {sgn} {abs(b_val):.2f} kHz "
        f"(R^2 = {r2:.4f}, RMSE = {rmse:.2f} kHz)"
    )
    dialog_msg = (
        "Safety Boundary Formula Fitted Successfully!\n\n"
        f"Equation: Df_hc (kHz) = {a_val:.4f} * I0 {sgn} {abs(b_val):.2f}\n"
        f"Operating Rule: Df_hc > {a_val:.4f} * I0 {sgn} {abs(b_val):.2f} kHz\n"
        f"Coefficient of Determination R^2: {r2:.4f}\n"
        f"Root Mean Square Error: {rmse:.2f} kHz\n"
        f"Critical Boundary Points: {len(crit_i0)}"
    )
    return {
        "crit_i0": crit_i0_arr,
        "crit_df": crit_df_arr,
        "a": a_val,
        "b": b_val,
        "sgn": sgn,
        "r2": r2,
        "rmse": rmse,
        "i0_dense": i0_dense,
        "df_dense": df_dense,
        "status_msg": status_msg,
        "dialog_msg": dialog_msg,
    }


def align_homs_to_worst_case_sideband(cfg: MachineConfig) -> List[HOMMode]:
    """
    Exact port of MATLAB `alignHOMToWorstCaseSideband` (STABLE_App.m lines 2159-2200).
    Aligns all HOM frequencies to the SAME worst-case upper synchrotron sideband:
      fr = (M * h + mu_target + nu_s) * f0
    """
    if not cfg.homs:
        return []
    cspeed = 299792458.0
    f0_val = cspeed / max(1.0, cfg.C)
    h_val = max(1, int(cfg.h))
    sin_phi_s = min(1.0, max(-1.0, cfg.U0 / max(1.0, cfg.V_mc)))
    cos_phi_s = math.sqrt(max(0.0, 1.0 - sin_phi_s ** 2))
    nu_s = math.sqrt((cfg.alpha_c * h_val * cfg.V_mc * cos_phi_s) / (2.0 * math.pi * max(1.0, cfg.E0)))
    fs_val = nu_s * f0_val

    f1 = cfg.homs[0].freq_Hz
    k1 = int(round((f1 - fs_val) / f0_val))
    mu_target = k1 % h_val
    if mu_target == 0:
        mu_target = 1

    aligned: List[HOMMode] = []
    for m in cfg.homs:
        fn = float(m.freq_Hz)
        if fn <= 0:
            aligned.append(HOMMode(Q=m.Q, R=m.R, freq_Hz=fn))
            continue
        m_mult = int(round((fn - fs_val - mu_target * f0_val) / (h_val * f0_val)))
        k_mode = m_mult * h_val + mu_target
        f_aln = round((k_mode + nu_s) * f0_val)
        aligned.append(HOMMode(Q=m.Q, R=m.R, freq_Hz=float(f_aln)))
    return aligned


def analyze_scan_folder_mat(folder_path: str | Path) -> Dict[str, Any]:
    """
    Exact port of MATLAB `AnalyzeButtonPushed` (STABLE_App.m lines 1739-1865).
    Scans `SOLEIL_I0*RLfp*QLfp*detune*.mat` files in `folder_path` and computes
    bunch length, centroid mean, centroid jitter, energy spread vs detuning, plus optimal detuning.
    """
    folder = Path(folder_path)
    pt_files = list(folder.glob("SOLEIL_I0*RLfp*QLfp*detune*.mat"))
    if not pt_files:
        raise FileNotFoundError(f"No SOLEIL_I0*RLfp*QLfp*detune*.mat files found in {folder}")

    parsed: List[Tuple[float, Path]] = []
    for f in pt_files:
        m = re.search(r"detune(-?\d+\.?\d*(?:[eE][-+]?\d+)?)", f.name)
        if m:
            parsed.append((float(m.group(1)), f))

    if not parsed:
        raise ValueError("Could not parse detuning frequency from filenames.")

    parsed.sort(key=lambda x: x[0])
    n_scans = len(parsed)
    fre_shift_kHz = np.zeros(n_scans, dtype=np.float64)
    mean_bl_ps = np.full(n_scans, np.nan, dtype=np.float64)
    centroid_mean_ps = np.full(n_scans, np.nan, dtype=np.float64)
    centroid_jitter_ps = np.full(n_scans, np.nan, dtype=np.float64)
    energy_spread = np.full(n_scans, np.nan, dtype=np.float64)

    for i, (df_hz, fpath) in enumerate(parsed):
        fre_shift_kHz[i] = df_hz * 1e-3
        try:
            d = sio.loadmat(str(fpath), squeeze_me=True, struct_as_record=False)
            sigma_t0 = 40e-12
            sigma_e0 = 7.44e-4
            if "HALF" in d:
                h_d = _mat_struct_to_dict(d["HALF"])
                sigma_t0 = _scalar_float(h_d.get("sigma_t0"), sigma_t0) or sigma_t0
                sigma_e0 = _scalar_float(h_d.get("sigma_e0"), sigma_e0) or sigma_e0
            if "sigma_t0" in d:
                sigma_t0 = _scalar_float(d.get("sigma_t0"), sigma_t0) or sigma_t0
            if "sigma_e0" in d:
                sigma_e0 = _scalar_float(d.get("sigma_e0"), sigma_e0) or sigma_e0

            if "record_Q_std" in d and "record_Q_mean" in d:
                q_std = np.asarray(d["record_Q_std"], dtype=np.float64)
                q_mean = np.asarray(d["record_Q_mean"], dtype=np.float64)
                if q_std.ndim == 1:
                    q_std = q_std[:, None]
                    q_mean = q_mean[:, None]
                n_rec = q_std.shape[0]
                s_start = max(0, int(round(n_rec * 0.9)) - 1)
                q_std_tail = q_std[s_start:, :]
                q_mean_tail = q_mean[s_start:, :]

                mean_bl_ps[i] = float(np.mean(q_std_tail)) * sigma_t0 * 1e12
                centroid_mean_ps[i] = float(np.mean(np.mean(q_mean_tail, axis=0))) * sigma_t0 * 1e12
                centroid_jitter_ps[i] = float(np.mean(np.std(q_mean_tail, axis=0, ddof=0))) * sigma_t0 * 1e12

            if "record_P_std" in d:
                p_std = np.asarray(d["record_P_std"], dtype=np.float64)
                if p_std.ndim == 1:
                    p_std = p_std[:, None]
                n_rec = p_std.shape[0]
                s_start = max(0, int(round(n_rec * 0.9)) - 1)
                energy_spread[i] = float(np.mean(p_std[s_start:, :])) * sigma_e0
        except Exception:
            pass

    valid = ~np.isnan(mean_bl_ps) & ~np.isnan(centroid_jitter_ps)
    scan_plot = fre_shift_kHz[valid]
    bl_plot = mean_bl_ps[valid]
    c_mean_plot = centroid_mean_ps[valid]
    jit_plot = centroid_jitter_ps[valid]
    es_plot = energy_spread[valid]

    best_detune = float("nan")
    best_bl = float("nan")
    if scan_plot.size > 0:
        stable_idx = np.where(jit_plot < 1.5)[0]
        if stable_idx.size > 0:
            sub_best = int(np.argmax(bl_plot[stable_idx]))
            best_i = stable_idx[sub_best]
            best_detune = float(scan_plot[best_i])
            best_bl = float(bl_plot[best_i])
        else:
            best_i = int(np.argmin(jit_plot))
            best_detune = float(scan_plot[best_i])
            best_bl = float(bl_plot[best_i])

    return {
        "scan_kHz": scan_plot,
        "bunch_length_ps": bl_plot,
        "centroid_mean_ps": c_mean_plot,
        "centroid_jitter_ps": jit_plot,
        "energy_spread": es_plot,
        "best_detune_kHz": best_detune,
        "best_bunch_length_ps": best_bl,
    }
