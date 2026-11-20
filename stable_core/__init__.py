"""
STABLE Core: Advanced Storage Ring RF Stability and Parameter Optimization Toolkit.

Provides:
- Instability classification (LLRF Nyquist, spatial DFT, Hilbert growth rate fitting, physical mode disambiguation)
- Physical stability penalty-constrained optimizer (RFOptimizer)
- Two-stage cascade accelerated parameter scanner (CascadeScanner)
- Safety boundary extractor and equation fitter (SafetyBoundaryExtractor)
- LLRF crescent stability window tuner (LLRFTuner)
"""

from stable_core.types import (
    StabilityStatus,
    InstabilityResult,
    MachineParameters,
)
from stable_core.classifier import InstabilityClassifier, calculate_d_mode_robinson
from stable_core.optimizer import RFOptimizer, OptimizationResult
from stable_core.cascade_scanner import CascadeScanner, ScanGridResult, ScanPointResult
from stable_core.boundary_fitter import SafetyBoundaryExtractor, BoundaryFitResult
from stable_core.llrf_tuner import LLRFTuner, LLRFTuningResult

__all__ = [
    "StabilityStatus",
    "InstabilityResult",
    "MachineParameters",
    "InstabilityClassifier",
    "calculate_d_mode_robinson",
    "RFOptimizer",
    "OptimizationResult",
    "CascadeScanner",
    "ScanGridResult",
    "ScanPointResult",
    "SafetyBoundaryExtractor",
    "BoundaryFitResult",
    "LLRFTuner",
    "LLRFTuningResult",
]

__version__ = "1.0.0"
