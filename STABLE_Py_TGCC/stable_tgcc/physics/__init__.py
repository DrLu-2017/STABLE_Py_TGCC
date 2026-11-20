"""
High-performance Multi-A100 GPU Physics & Longitudinal Tracking Engine.
Faithfully ports MATLAB STABLE (machine.m, MC_Cavity_Control.m, DFB_Set.m, MPMBTrackBunchLengthingwoHighQ.m).
"""

from stable_tgcc.physics.backend import ArrayBackend, detect_available_gpus
from stable_tgcc.physics.machine_init import (
    build_half_machine,
    calc_bunch_form_factor_eq,
    detune_hc_calc,
    compute_optimal_mc_detune,
    align_hom_to_worst_sideband,
    generate_fill_pattern,
)
from stable_tgcc.physics.cavity_feedback import (
    dfb_set,
    dfb_calc_properties,
    mc_cavity_control,
    voltage_load_calc,
)
from stable_tgcc.physics.analysis import (
    calculate_d_mode_robinson,
    analyze_llrf_nyquist,
    classify_tracking_result,
)
from stable_tgcc.physics.tracking_engine import run_single_tracking
from stable_tgcc.physics.multi_gpu_scanner import run_2d_stability_scan

__all__ = [
    "ArrayBackend",
    "detect_available_gpus",
    "build_half_machine",
    "calc_bunch_form_factor_eq",
    "detune_hc_calc",
    "compute_optimal_mc_detune",
    "align_hom_to_worst_sideband",
    "generate_fill_pattern",
    "dfb_set",
    "dfb_calc_properties",
    "mc_cavity_control",
    "voltage_load_calc",
    "calculate_d_mode_robinson",
    "analyze_llrf_nyquist",
    "classify_tracking_result",
    "run_single_tracking",
    "run_2d_stability_scan",
]
