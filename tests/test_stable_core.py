"""
Unit test suite for STABLE core modules.
"""

import unittest
import numpy as np

from stable_core.types import StabilityStatus, MachineParameters, InstabilityResult
from stable_core.classifier import InstabilityClassifier, calculate_d_mode_robinson
from stable_core.optimizer import RFOptimizer
from stable_core.cascade_scanner import CascadeScanner
from stable_core.boundary_fitter import SafetyBoundaryExtractor
from stable_core.llrf_tuner import LLRFTuner


class TestStableCore(unittest.TestCase):
    """Test suite for the four refactored modules and core classes."""

    def setUp(self) -> None:
        self.machine = MachineParameters()
        self.classifier = InstabilityClassifier(self.machine)

    def test_types_and_defaults(self) -> None:
        """Verify SOLEIL II worst-case RF (IDs open) machine parameter calculations."""
        self.assertAlmostEqual(self.machine.C, 354.097, places=3)
        self.assertAlmostEqual(self.machine.Q_mc, 6000.0, places=1)
        self.assertAlmostEqual(self.machine.tau_z, 13.7e-3, places=4)
        self.assertAlmostEqual(self.machine.sigma_delta0, 0.92e-3, places=5)
        self.assertAlmostEqual(self.machine.f_rf, 352.2e6, delta=5e5)
        self.assertGreater(self.machine.T0, 0.0)
        self.assertAlmostEqual(self.machine.T0 * self.machine.c, self.machine.C, delta=1.0)

    def test_classifier_nyquist(self) -> None:
        """Verify physical frequency-domain LLRF Nyquist pole and margin analysis (including SOLEIL II 2.0 us delay)."""
        # SOLEIL II actual operating delay: 704 buckets (~2.0 us = 1.69 * T0) at G_dfb = 0.2
        is_s2_stable, s2_margin, _ = self.classifier.analyze_llrf_nyquist(
            K_p=0.2, delay_sec=self.machine.dfb_delay
        )
        self.assertTrue(is_s2_stable)
        self.assertGreater(s2_margin, 0.5)

        # Standard half-turn loop delay tau_d ~ T0 / 2 at K_p = 1.0
        delay_stable = self.machine.T0 / 2.0
        is_stable, margin, dist = self.classifier.analyze_llrf_nyquist(K_p=1.0, delay_sec=delay_stable)
        self.assertTrue(is_stable)
        self.assertGreater(margin, 0.05)

        # Excessive gain & delay exceeds phase-crossover G_max -> encirclement instability
        is_unstable, margin_bad, _ = self.classifier.analyze_llrf_nyquist(K_p=5.0, delay_sec=1.9 * self.machine.T0)
        self.assertFalse(is_unstable)
        self.assertEqual(margin_bad, 0.0)

    def test_iterative_form_factor_toggle(self) -> None:
        """Verify self-consistent 4th-harmonic bunch form factor F_4(sigma_tau) iteration & toggle."""
        df_point = self.machine.get_nominal_hc_detuning(500.0, target_xi=1.0, use_form_factor=False)
        df_iter = self.machine.get_nominal_hc_detuning(500.0, target_xi=1.0, use_form_factor=True)
        # Stretched 4th-harmonic form factor F_4 < 1 reduces V_hc, requiring ~6-7% lower detuning for xi_eff = 1.0
        self.assertLess(df_iter, df_point)
        self.assertGreater(df_point - df_iter, 4.0)

        eq_on = self.classifier.solve_equilibrium_form_factor(0.5, 88e3, use_form_factor=True)
        eq_off = self.classifier.solve_equilibrium_form_factor(0.5, 88e3, use_form_factor=False)
        self.assertEqual(eq_off["F_hc"], 1.0)
        self.assertLess(eq_on["F_hc"], 0.96)
        self.assertGreater(eq_on["F_hc"], 0.85)
        self.assertGreater(eq_on["sigma_t_eq"], 3.5 * eq_on["sigma_t0"])

    def test_classifier_semianalytical(self) -> None:
        """Verify semi-analytical evaluation rules."""
        # Under nominal beam loading and detuning (xi ~ 0.90)
        res_stable = self.classifier.evaluate_semianalytical(I0_mA=200.0, df_khz=40.0)
        self.assertEqual(res_stable.status, StabilityStatus.STABLE)
        self.assertTrue(res_stable.is_stable)

        # High current and over-stretching triggers PTBL (I0 >= 365 mA, xi >= 1.05)
        res_ptbl = self.classifier.evaluate_semianalytical(I0_mA=500.0, df_khz=50.0)
        self.assertEqual(res_ptbl.status, StabilityStatus.PTBL)
        self.assertFalse(res_ptbl.is_stable)

    def test_optimizer_stability_constraint(self) -> None:
        """Verify optimizer finds stable operating point under penalty constraint."""
        opt = RFOptimizer(self.machine, self.classifier)
        res = opt.optimize_harmonic_cavity_detuning(
            I0_mA=400.0,
            initial_detuning_khz=80.0,
            detuning_bounds_khz=(40.0, 120.0),
            method="COBYLA",
        )
        self.assertTrue(res.success)
        self.assertEqual(res.status, StabilityStatus.STABLE)
        self.assertGreater(res.bunch_length_factor, 1.5)
        self.assertGreaterEqual(res.xi, 0.70)
        self.assertLessEqual(res.xi, 1.05)
        self.assertIn("detuning_khz", res.optimal_params)

    def test_cascade_scanner(self) -> None:
        """Verify two-stage cascade pruning and parallel execution."""
        scanner = CascadeScanner(self.machine, max_workers=2)
        x_vals = np.linspace(30.0, 60.0, 4)
        y_vals = np.linspace(200.0, 500.0, 4)

        grid = scanner.scan_2d(
            x_name="detuning_khz",
            x_values=x_vals,
            y_name="I0_mA",
            y_values=y_vals,
            fixed_params={"K_p": 1.0},
        )
        self.assertEqual(grid.status_grid.shape, (4, 4))
        self.assertGreaterEqual(grid.prune_ratio, 0.0)
        self.assertGreaterEqual(grid.estimated_speedup, 1.0)

    def test_boundary_fitter(self) -> None:
        """Verify boundary extractor and equation fitting."""
        extractor = SafetyBoundaryExtractor()

        # Synthetic threshold: y_crit(x) = 0.1 * x + 5.0
        def eval_func(x: float, y: float) -> bool:
            return y > (0.1 * x + 5.0)

        fit = extractor.extract_from_evaluator(
            eval_func=eval_func,
            x_name="I0_mA",
            x_values=[100.0, 200.0, 300.0, 400.0, 500.0],
            y_name="detuning_khz",
            y_bounds=(0.0, 80.0),
            degree=1,
            safe_region="above",
        )
        self.assertAlmostEqual(fit.coefficients[0], 0.1, places=2)
        self.assertAlmostEqual(fit.coefficients[1], 5.0, places=1)
        self.assertGreater(fit.r2_score, 0.99)
        self.assertIn("detuning_khz", fit.inequality_str)
        self.assertTrue(fit.is_safe(300.0, 40.0))
        self.assertFalse(fit.is_safe(300.0, 20.0))

    def test_llrf_tuner(self) -> None:
        """Verify LLRF crescent stability window and Chebyshev center."""
        tuner = LLRFTuner(self.machine, self.classifier)
        res = tuner.scan_feedback_plane(
            kp_range=(0.2, 3.0),
            n_kp=15,
            n_delay=15,
            I0_mA=400.0,
            df_khz=85.0,
        )
        self.assertGreater(res.stable_area_fraction, 0.0)
        self.assertGreater(res.optimal_kp, 0.0)
        self.assertGreater(res.optimal_delay, 0.0)
        self.assertGreater(res.optimal_margin_score, 0.0)

    def test_mode0_vs_lcbi_classification(self) -> None:
        """Verify pure Mode-0 (mu=0) oscillation is never misclassified as LCBI (mu != 0)."""
        n_turns, h = 500, 416
        t = np.arange(n_turns)[:, None]
        bunches = np.arange(h)[None, :]

        # Case 1: Large pure Mode-0 (mu=0) oscillation near flat potential (350 mA, 70 kHz, xi ~ 0.99)
        q_mode0 = 8.0 * np.sin(2.0 * np.pi * 0.02 * t) * np.ones((1, h))
        p_mode0 = 0.8 * np.cos(2.0 * np.pi * 0.02 * t) * np.ones((1, h))
        q_std = np.full((n_turns, h), 4.0)
        p_std = np.full((n_turns, h), 1.45)

        res_m0 = self.classifier.evaluate_from_signals(
            q_mode0, p_mode0, q_std, p_std, I0=0.35, fre_shift_hc=70e3
        )
        self.assertNotEqual(res_m0.status, StabilityStatus.LCBI)
        self.assertEqual(res_m0.status, StabilityStatus.COUPLED_DQ)
        self.assertEqual(res_m0.mode_id, 0)

        # Case 2: Genuine coupled-bunch mode mu=71 driven by HOM
        phase_cb = 2.0 * np.pi * (0.02 * t + 71.0 * bunches / h)
        q_lcbi = 4.0 * np.exp(0.003 * t) * np.sin(phase_cb)
        p_lcbi = 0.5 * np.exp(0.003 * t) * np.cos(phase_cb)

        res_cb = self.classifier.evaluate_from_signals(
            q_lcbi, p_lcbi, q_std, p_std, I0=0.35, fre_shift_hc=90e3
        )
        self.assertEqual(res_cb.status, StabilityStatus.LCBI)
        self.assertIn(res_cb.mode_id, (71, h - 71))

    def test_calculate_d_mode_robinson_half_pshc(self) -> None:
        """Verify analytical D-mode Robinson frequency, growth rate, and threshold on HALF PSHC."""
        # HALF storage ring superconducting PSHC benchmark (T. He et al., PRAB 2024)
        dm_unstable = calculate_d_mode_robinson(
            I0=0.35, E_GeV=2.2, C_m=479.86, alpha_c=8.1e-5, tau_z_s=13.5e-3,
            U0_keV=380.0, Vrf_MV=1.25, R1_Q1=90.0, Q1=2e5, df_r1_Hz=-12e3,
            n_harmonic=3, Rn_Qn=44.5, Qn=2e8, df_rn_Hz=8.0e3, h=800,
        )
        self.assertTrue(dm_unstable["is_unstable"])
        self.assertGreater(dm_unstable["growth_rate_s1"], 50.0)
        self.assertAlmostEqual(
            dm_unstable["f_Dmode_Hz"] + dm_unstable["df1_Hz"], 8000.0, places=2
        )

        dm_stable = calculate_d_mode_robinson(
            I0=0.35, E_GeV=2.2, C_m=479.86, alpha_c=8.1e-5, tau_z_s=13.5e-3,
            U0_keV=380.0, Vrf_MV=1.25, R1_Q1=90.0, Q1=2e5, df_r1_Hz=-12e3,
            n_harmonic=3, Rn_Qn=44.5, Qn=2e8, df_rn_Hz=30e3, h=800,
        )
        self.assertFalse(dm_stable["is_unstable"])
        self.assertLess(dm_stable["growth_rate_s1"], 0.0)
        self.assertGreater(dm_stable["thresh_detune_Hz"], 7.0e3)
        self.assertLess(dm_stable["thresh_detune_Hz"], 16.5e3)
        self.assertGreater(dm_stable["opt_detune_Hz"], dm_stable["thresh_detune_Hz"])
        self.assertTrue(dm_stable["is_stretching_safe"])

    def test_dmode_threshold_curve_and_cascade_buffer(self) -> None:
        """Verify threshold curve scan and Stage-1 critical-zone buffer in CascadeScanner."""
        half_machine = MachineParameters(
            C=479.86, E0=2.2e9, U0=380e3, h=800, alpha_c=8.1e-5, tau_z=13.5e-3,
            V_mc=1.25e6, Q_mc=2e4, R_mc=1.8e6, RoQ_mc=90.0, fre_shift_mc=-50e3,
            n_hc=3, Q_hc=2e8, R_hc=8.9e9, RoQ_hc=44.5, fre_shift_hc=30e3,
            DFB_enable=False,
        )
        clf_half = InstabilityClassifier(half_machine)
        curve = clf_half.compute_dmode_threshold_curve(np.array([150.0, 250.0, 350.0]))
        self.assertEqual(len(curve["thresh_khz"]), 3)
        self.assertTrue(np.all(np.diff(curve["thresh_khz"]) > 0.0))
        self.assertTrue(np.all(curve["opt_khz"] > curve["thresh_khz"]))

        # Stage-1 cascade check at low current (I0 = 10 mA) where D-mode threshold (~5.0 kHz)
        # exceeds flat-potential optimal detuning (~3.4 kHz), so xi < 0.85 (pure single-well D-mode!):
        # deep D-mode (0.42 * thresh) is pruned in Stage 1, while near-threshold (0.99 * thresh) passes to Stage 2!
        half_highq = MachineParameters(
            C=479.86, E0=2.2e9, U0=380e3, h=800, alpha_c=8.1e-5, tau_z=13.5e-3,
            V_mc=1.25e6, Q_mc=2e4, R_mc=1.8e6, RoQ_mc=90.0, fre_shift_mc=-50e3,
            n_hc=3, Q_hc=1e9, R_hc=8.9e10, RoQ_hc=89.0, fre_shift_hc=30e3,
            DFB_enable=False,
        )
        clf_hq = InstabilityClassifier(half_highq)
        dm_10mA = clf_hq.evaluate_dmode_robinson(0.010, 15e3, compute_threshold=True)
        self.assertFalse(dm_10mA["is_stretching_safe"])  # At 10 mA, thresh (~5.0 kHz) > opt (~3.4 kHz)
        thr_10mA_khz = dm_10mA["thresh_detune_Hz"] * 1e-3
        scanner = CascadeScanner(half_highq, max_workers=1)
        grid = scanner.scan_2d(
            x_name="detuning_khz",
            x_values=np.array([0.42 * thr_10mA_khz, 0.99 * thr_10mA_khz, 20.0]),
            y_name="I0_mA",
            y_values=np.array([10.0]),
            fixed_params={"K_p": 0.0, "critical_margin_s1": 60.0},
        )
        self.assertEqual(grid.stage_grid[0, 0], 1)  # Manifestly unstable D-mode -> pruned at Stage 1
        self.assertEqual(grid.stage_grid[0, 1], 2)  # Near-threshold (0 < Omega_i < crit_margin) -> passed to Stage 2
        self.assertEqual(grid.stage_grid[0, 2], 2)  # Stable (20 kHz) -> passed to Stage 2

    def test_mode0_spectrum_benchmark(self) -> None:
        """Verify Mode-0 coherent FFT peak matches analytical D-mode frequency Omega_r / 2pi."""
        half_machine = MachineParameters(
            C=479.86, E0=2.2e9, U0=380e3, h=800, alpha_c=8.1e-5, tau_z=13.5e-3,
            V_mc=1.25e6, Q_mc=2e5, R_mc=1.8e7, RoQ_mc=90.0, fre_shift_mc=-12e3,
            n_hc=3, Q_hc=2e8, R_hc=8.9e9, RoQ_hc=44.5, fre_shift_hc=15e3,
            DFB_enable=False,
        )
        clf_half = InstabilityClassifier(half_machine)
        dm = clf_half.evaluate_dmode_robinson(0.35, 15e3, compute_threshold=False)
        f_theory = dm["f_Dmode_Hz"]

        n_turns, n_bun = 2000, 32
        t_sec = np.arange(n_turns)[:, None] * half_machine.T0
        q_sig = 3.0 * np.cos(2.0 * np.pi * f_theory * t_sec) * np.ones((1, n_bun))
        p_sig = -0.3 * np.sin(2.0 * np.pi * f_theory * t_sec) * np.ones((1, n_bun))
        q_std = np.full((n_turns, n_bun), 3.2)
        p_std = np.full((n_turns, n_bun), 1.0)

        res = clf_half.evaluate_from_signals(
            q_sig, p_sig, q_std, p_std, I0=0.35, fre_shift_hc=15e3, recor_step=1
        )
        self.assertAlmostEqual(res.details["f_mode0_Hz"], f_theory, delta=120.0)


if __name__ == "__main__":
    unittest.main()
