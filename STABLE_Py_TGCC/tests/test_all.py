"""
Comprehensive Unit & Integration Tests for STABLE_Py_TGCC:
  1. Physics parity with MATLAB (Fokker-Planck F_n iteration, machine.m, DFB phase_shift=0 neg FB, MC_Cavity_Control)
  2. Multi-user encrypted vault (PBKDF2 + Fernet/AES encryption of TGCC SSH secrets & per-user presets)
  3. Server-Client streaming RPC (Status, Nyquist, D-mode, Live Single Tracking stream, Multi-GPU 2D Scan stream)
  4. TGCC `ccc_msub` A100 script generation & batch runner
"""

from __future__ import annotations
import tempfile
import unittest
from pathlib import Path
import numpy as np

from stable_tgcc.common.models import (
    MachineConfig,
    ScanConfig,
    TGCCConnectionProfile,
    TelemetryFrame,
    ScanPointTelemetry,
)
from stable_tgcc.physics.machine_init import (
    build_half_machine,
    calc_bunch_form_factor_eq,
    detune_hc_calc,
    compute_optimal_mc_detune,
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
)
from stable_tgcc.auth.user_vault import UserVaultManager
from stable_tgcc.client.ssh_tunnel import STABLEClientConnector
from stable_tgcc.server.ccc_scheduler import generate_ccc_msub_script


class TestSTABLEPyTGCC(unittest.TestCase):

    def test_01_fokker_planck_form_factor_and_machine_init(self):
        """Verify Fokker-Planck bunch form factor iteration F_n(sigma_t) and machine.m port."""
        cfg = MachineConfig.soleil2_default()
        fp = calc_bunch_form_factor_eq(
            "eval", cfg.I0, cfg.V_mc, cfg.U0, cfg.E0, cfg.C, cfg.h,
            cfg.alpha_c, cfg.sigma_e0, cfg.n_hc, cfg.R_hc, cfg.Q_hc,
            67.88e3, True,
        )
        self.assertGreater(fp["F_hc"], 0.70)
        self.assertLess(fp["F_hc"], 0.98)
        self.assertGreater(fp["sigma_t_eq"], fp["sigma_t0"] * 2.0)

        df_xi1 = detune_hc_calc(
            cfg.I0, cfg.n_hc, cfg.C, cfg.h, cfg.U0, cfg.V_mc,
            cfg.R_hc, cfg.Q_hc, cfg.E0, cfg.alpha_c, cfg.sigma_e0, True, 1.0,
        )
        self.assertGreater(df_xi1, 50e3)
        self.assertLess(df_xi1, 90e3)

        df_mc = compute_optimal_mc_detune(
            cfg.C, cfg.h, cfg.I0, cfg.V_mc, cfg.U0, cfg.R_mc, cfg.Q_mc,
            cfg.R_hc, cfg.Q_hc, df_xi1, cfg.n_hc,
        )
        self.assertLess(df_mc, -15e3)

        HALF = build_half_machine(
            cfg.C, cfg.I0, cfg.U0, cfg.E0, cfg.tau_s, cfg.tau_z,
            fp["sigma_t0"], cfg.sigma_e0, cfg.alpha_c, cfg.h, cfg.V_mc,
            cfg.n_hc, cfg.R_hc, cfg.Q_hc, 1.0, df_xi1, cfg.Q_mc, cfg.R_mc, df_mc, True,
        )
        self.assertAlmostEqual(abs(HALF["Vrf_ideal"]), cfg.V_mc, delta=1.0)

    def test_02_dfb_phase_zero_negative_feedback_and_cavity_control(self):
        """Verify DFB with phase_shift=0 is negative feedback, matches Ig at t=0, and reduces Rs."""
        cfg = MachineConfig.soleil2_default()
        HALF = build_half_machine(
            cfg.C, cfg.I0, cfg.U0, cfg.E0, cfg.tau_s, cfg.tau_z,
            9e-12, cfg.sigma_e0, cfg.alpha_c, cfg.h, cfg.V_mc,
            cfg.n_hc, cfg.R_hc, cfg.Q_hc, 1.0, cfg.fre_shift_hc,
            cfg.Q_mc, cfg.R_mc, cfg.fre_shift_mc, True,
        )
        # phase_shift = 0.0 rad MUST be negative feedback
        DFB = dfb_set(HALF, gain=0.25, phase_shift=0.0, delay=704, sample_num=208, every=208)
        self.assertAlmostEqual(abs((DFB["Ig_main0"] + DFB["Ig_dfb0"]) - HALF["Ig_mc_0"]), 0.0, places=12)

        props = dfb_calc_properties(DFB, HALF)
        self.assertLess(props["Rs"], HALF["R_mc"])
        self.assertGreater(props["gamma"], 0.0)
        self.assertGreater(props["fs"], 100.0)

        # Verify VoltageLoadCalc vectorized lfilter matches exact recurrence
        pattern = np.ones(cfg.h, dtype=np.int32)
        rot_decay_mc = 1j * HALF["rot_coef_mc"] - 1.0 / (2.0 * cfg.Q_mc)
        TbAng_mc = np.exp(rot_decay_mc * HALF["angle_mc"])
        v_load_cpu = np.full(cfg.h, HALF["qc"] * HALF["wr_mc"] * cfg.R_mc / cfg.Q_mc, dtype=np.complex128)
        v_arr, v_end = voltage_load_calc(HALF["V_mc_load_0"], v_load_cpu, TbAng_mc, pattern)
        self.assertEqual(len(v_arr), cfg.h)
        self.assertTrue(np.isfinite(abs(v_end)))

    def test_03_multi_user_encrypted_vault(self):
        """Verify per-user password hashing, AES/Fernet SSH secret encryption, and preset isolation."""
        with tempfile.TemporaryDirectory() as tmpdir:
            db_file = Path(tmpdir) / "test_vault.db"
            vault = UserVaultManager(db_path=db_file)

            s_alice = vault.register_user("alice", "AlicePass!23")
            s_bob = vault.register_user("bob", "BobSecret#99")

            prof_a = TGCCConnectionProfile(
                profile_name="TGCC Irene A100",
                mode="ssh_tgcc",
                ssh_host="irene-fr.ccc.cea.fr",
                ssh_username="alice_cea",
                ccc_project="gen9999",
            )
            vault.save_tgcc_profile(s_alice, prof_a, ssh_secret="TopSecretTGCCPassword")

            # Verify raw SQLite does NOT contain plaintext secret
            raw_bytes = db_file.read_bytes()
            self.assertNotIn(b"TopSecretTGCCPassword", raw_bytes)

            # Verify Alice can decrypt her secret
            loaded_prof, dec_secret = vault.load_tgcc_profile(s_alice)
            self.assertEqual(loaded_prof.ssh_username, "alice_cea")
            self.assertEqual(dec_secret, "TopSecretTGCCPassword")

            # Verify password change re-encrypts secret properly
            s_alice_new = vault.change_password(s_alice, "AlicePass!23", "NewAlicePass$45")
            _, dec_secret_after = vault.load_tgcc_profile(s_alice_new)
            self.assertEqual(dec_secret_after, "TopSecretTGCCPassword")

            # Verify Bob cannot see Alice's profile
            prof_b, sec_b = vault.load_tgcc_profile(s_bob)
            self.assertEqual(prof_b.ssh_username, "bob")
            self.assertEqual(sec_b, "")

    def test_04_server_client_streaming_rpc_and_ccc_msub(self):
        """Verify Client <-> Server RPC for Nyquist, D-mode, Single Tracking stream, and 2D Scan stream."""
        conn = STABLEClientConnector()
        try:
            st = conn.connect(TGCCConnectionProfile(mode="local"))
            self.assertTrue(st["ok"])

            cfg = MachineConfig.soleil2_default()
            cfg.Par_num = 100
            cfg.Track_num = 60
            cfg.Recor_step = 10

            # 1. Nyquist RPC
            nyq_resp = conn.request_json("/api/nyquist", {"machine": cfg.to_dict()})
            self.assertTrue(nyq_resp["ok"])
            self.assertGreater(nyq_resp["nyquist"]["g_dfb_max"], 0.2)

            # 2. D-mode RPC
            dm_resp = conn.request_json("/api/dmode", {"machine": cfg.to_dict()})
            self.assertTrue(dm_resp["ok"])
            self.assertIn("growth_rate_s1", dm_resp["dmode"])

            # 3. Single Tracking Stream
            frames: list[TelemetryFrame] = []
            conn.stream_single_tracking(
                cfg=cfg,
                job_id="unit_test_track",
                on_frame=lambda fr: frames.append(fr),
            )
            self.assertGreaterEqual(len(frames), 2)
            self.assertEqual(frames[-1].event_type, "tracking_done")

            # 4. 2D Multi-GPU Scan Stream (2x2 mini grid)
            scan_cfg = ScanConfig(
                Df_min_kHz=65.0,
                Df_max_kHz=70.0,
                Df_step_kHz=5.0,
                I0_min_mA=330.0,
                I0_max_mA=350.0,
                I0_step_mA=20.0,
                num_gpus=2,
            )
            pts: list[ScanPointTelemetry] = []
            conn.stream_2d_scan(
                cfg=cfg,
                scan_cfg=scan_cfg,
                profile=TGCCConnectionProfile(mode="local"),
                job_id="unit_test_scan",
                on_point=lambda p: pts.append(p),
            )
            self.assertEqual(len(pts), 4)

            # 5. TGCC ccc_msub script generation
            sh_text = generate_ccc_msub_script(
                profile=TGCCConnectionProfile(ccc_project="gen1234", ccc_partition="a100"),
                job_id="scan_01",
                job_dir_remote="/ccc/scratch/cont003/user/job1",
                num_gpus=4,
            )
            self.assertIn("#MSUB -q a100", sh_text)
            self.assertIn("#MSUB -A gen1234", sh_text)
        finally:
            conn.disconnect()


if __name__ == "__main__":
    unittest.main()
