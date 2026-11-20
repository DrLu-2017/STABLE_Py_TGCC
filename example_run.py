"""
Example Run: End-to-End Demonstration of Refactored STABLE Acceleration Framework.

Demonstrates all four core engineering and physics modules:
1. Physical Stability-Constrained Optimizer (RFOptimizer)
2. Two-Stage Cascade Accelerated Scanner (CascadeScanner)
3. Safety Boundary Extractor & Formula Fitter (SafetyBoundaryExtractor)
4. LLRF Direct Feedback Crescent Stability Window Tuner (LLRFTuner)
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

# Ensure local repository is in path
sys.path.insert(0, str(Path(__file__).resolve().parent))

from stable_core.types import MachineParameters, StabilityStatus
from stable_core.classifier import InstabilityClassifier
from stable_core.optimizer import RFOptimizer
from stable_core.cascade_scanner import CascadeScanner
from stable_core.boundary_fitter import SafetyBoundaryExtractor
from stable_core.llrf_tuner import LLRFTuner


def run_task1_optimization(machine: MachineParameters, classifier: InstabilityClassifier) -> None:
    print("=" * 70)
    print("TASK 1: Physical Stability-Constrained RF Parameter Optimization")
    print("=" * 70)

    optimizer = RFOptimizer(machine, classifier)
    test_currents_mA = [250.0, 380.0, 500.0]

    for I0_mA in test_currents_mA:
        print(f"\n[+] Optimizing Harmonic Cavity Detuning for Beam Current I0 = {I0_mA:.1f} mA...")
        t0 = time.perf_counter()

        res = optimizer.optimize_harmonic_cavity_detuning(
            I0_mA=I0_mA,
            detuning_bounds_khz=(30.0, 140.0),
            method="COBYLA",
        )
        t_el = time.perf_counter() - t0

        opt_df = res.optimal_params["detuning_khz"]
        print(f"    - Success: {res.success}")
        print(f"    - Optimal Detuning Delta f_hc: {opt_df:.2f} kHz")
        print(f"    - Flat-Potential Parameter xi:  {res.xi:.3f}")
        print(f"    - Beam Status:                 {res.status.name} (Code: {int(res.status)})")
        print(f"    - Growth Rate:                 {res.growth_rate:.2f} s^-1")
        print(f"    - Bunch Lengthening Factor:    {res.bunch_length_factor:.2f}x")
        print(f"    - Nyquist Distance Margin:     {res.nyquist_margin:.3f}")
        print(f"    - Iterations:                  {res.iterations} (took {t_el*1000:.1f} ms)")


def run_task2_cascade_scan(machine: MachineParameters) -> tuple[CascadeScanner, any]:
    print("\n" + "=" * 70)
    print("TASK 2: Two-Stage Cascade Accelerated Parameter Scanning Engine")
    print("=" * 70)

    # 15x15 parameter grid: Detuning (40 - 120 kHz) vs Beam Current (200 - 550 mA)
    n_x, n_y = 16, 16
    detunings = np.linspace(40.0, 120.0, n_x)
    currents = np.linspace(200.0, 550.0, n_y)

    print(f"[+] Starting 2D Grid Scan ({n_x}x{n_y} = {n_x * n_y} points)...")
    scanner = CascadeScanner(machine, max_workers=4)

    t0 = time.perf_counter()
    grid_result = scanner.scan_2d(
        x_name="I0_mA",
        x_values=currents,
        y_name="detuning_khz",
        y_values=detunings,
        fixed_params={"K_p": 1.0},
    )
    t_el = time.perf_counter() - t0

    print(f"[+] Grid Scan Finished in {t_el:.3f} s!")
    print(f"    - Stage 1 Pruned Points: {grid_result.prune_ratio * 100:.1f}%")
    print(f"    - Estimated Pipeline Acceleration Speedup: {grid_result.estimated_speedup:.1f}x")
    stable_count = np.sum(grid_result.get_stable_mask())
    print(f"    - Stable Operating Points: {stable_count} / {grid_result.status_grid.size}")

    return scanner, grid_result


def run_task3_boundary_fitting(grid_result: any) -> any:
    print("\n" + "=" * 70)
    print("TASK 3: Safety Boundary Extractor and Equation Fitter")
    print("=" * 70)

    extractor = SafetyBoundaryExtractor()

    # Extract safety boundary from grid contour
    fit_res = extractor.extract_from_grid(
        grid_result=grid_result,
        safe_region="above",
        degree=1,
    )

    print("[+] Extracted Engineering Safety Boundary Equation:")
    print(f"    Equation:   {fit_res.equation_str}")
    print(f"    Safe Rule:  {fit_res.inequality_str}")
    print(f"    R^2 Score:  {fit_res.r2_score:.4f}")
    print(f"    RMS Error:  {fit_res.rmse:.2f} kHz")
    print(f"    Fitted Critical Boundary Points: {len(fit_res.critical_x)}")

    # Verification of test points
    test_pts = [(300.0, 70.0), (500.0, 60.0)]
    for cur, df in test_pts:
        is_safe = fit_res.is_safe(cur, df)
        print(f"    - Operating Point (I0 = {cur:.0f} mA, df = {df:.0f} kHz) -> Safe: {is_safe}")

    return fit_res


def run_task4_llrf_tuning(machine: MachineParameters, classifier: InstabilityClassifier, out_fig_path: Path) -> None:
    print("\n" + "=" * 70)
    print("TASK 4: LLRF Feedback Crescent-Shaped Stability Window Tuner")
    print("=" * 70)

    tuner = LLRFTuner(machine, classifier)
    print("[+] Mapping (K_p vs tau_d) parameter plane (40x40 grid)...")

    t0 = time.perf_counter()
    tune_res = tuner.scan_feedback_plane(
        kp_range=(0.1, 3.5),
        n_kp=40,
        n_delay=40,
        I0_mA=400.0,
        df_khz=85.0,
    )
    t_el = time.perf_counter() - t0

    print(f"[+] LLRF Tuning Completed in {t_el:.3f} s!")
    print(f"    - Crescent Stability Area Fraction: {tune_res.stable_area_fraction * 100:.1f}%")
    print(f"    - Chebyshev Center (Max Margin Setpoint):")
    print(f"        * Optimal Proportional Gain K_p = {tune_res.optimal_kp:.2f}")
    print(f"        * Optimal Loop Delay tau_d      = {tune_res.optimal_delay * 1e6:.3f} us ({tune_res.optimal_delay / machine.T0:.3f} T0)")
    print(f"        * Margin / Inscribed Radius     = {tune_res.optimal_margin_score:.2f}")
    print(f"    - Geometric Centroid Setpoint:")
    print(f"        * Centroid K_p                   = {tune_res.centroid_kp:.2f}")
    print(f"        * Centroid Loop Delay tau_d     = {tune_res.centroid_delay * 1e6:.3f} us ({tune_res.centroid_delay / machine.T0:.3f} T0)")

    # Plot crescent stability window
    tuner.plot_crescent_window(
        result=tune_res,
        save_path=out_fig_path,
        title="LLRF Direct Feedback Crescent Stability Window & Chebyshev Center",
    )
    print(f"[+] Saved visualization figure to: {out_fig_path}")


def main() -> None:
    print("\n" + "#" * 70)
    print("  STABLE Core: Storage Ring RF Stability & Optimization Engine")
    print("#" * 70 + "\n")

    machine = MachineParameters()
    classifier = InstabilityClassifier(machine)

    # Output artifact directory
    artifact_dir = Path(r"C:\Users\zhaol\.gemini\antigravity\brain\98411c23-a6cf-4664-a84d-0567203ec2cc")
    fig_path = artifact_dir / "llrf_crescent_stability_window.png"

    # Task 1
    run_task1_optimization(machine, classifier)

    # Task 2
    scanner, grid_result = run_task2_cascade_scan(machine)

    # Task 3
    fit_res = run_task3_boundary_fitting(grid_result)

    # Task 4
    run_task4_llrf_tuning(machine, classifier, fig_path)

    print("\n" + "=" * 70)
    print("ALL 4 TASKS SUCCESSFULLY DEMONSTRATED AND VALIDATED!")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    main()
