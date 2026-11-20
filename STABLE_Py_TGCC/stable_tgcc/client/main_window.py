"""
Full PySide6 Desktop GUI Client for STABLE_Py_TGCC.
Replicates all parameter tabs, analytical diagnostics, and the exact MATLAB `plotStabilityMap`
visualization of STABLE_App.m while streaming live per-turn & per-point data from TGCC Topaze A100.
"""

from __future__ import annotations
from typing import Dict, Any, List, Optional, Tuple, Callable
import math
import time
import threading
from pathlib import Path
import numpy as np
import scipy.io as sio
from scipy.interpolate import RectBivariateSpline
from scipy.ndimage import convolve

from PySide6.QtCore import Qt, Signal, QObject, QTimer, QSize
from PySide6.QtGui import QFont, QGuiApplication
from PySide6.QtWidgets import (
    QMainWindow,
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QFormLayout,
    QSplitter,
    QTabWidget,
    QGroupBox,
    QLabel,
    QLineEdit,
    QDoubleSpinBox,
    QSpinBox,
    QCheckBox,
    QComboBox,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QTextEdit,
    QProgressBar,
    QMessageBox,
    QInputDialog,
    QFileDialog,
    QHeaderView,
    QSizePolicy,
    QScrollArea,
    QFrame,
)

from matplotlib.backends.backend_qtagg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.figure import Figure

from stable_tgcc.auth.user_vault import UserVaultManager, UserSession
from stable_tgcc.common.models import (
    HOMMode,
    MachineConfig,
    ScanConfig,
    TGCCConnectionProfile,
    TelemetryFrame,
    ScanPointTelemetry,
)
from stable_tgcc.client.busy_overlay import AsyncTaskProgressDialog, SpinningRingWidget
from stable_tgcc.client.matlab_compat import (
    restore_scan_config_from_mat,
    load_config_from_mat_or_m,
    machine_config_to_matlab_struct,
    find_latest_default_scan_summary,
    fit_safety_boundary,
    align_homs_to_worst_case_sideband,
    analyze_scan_folder_mat,
)
from stable_tgcc.client.ssh_tunnel import STABLEClientConnector
from stable_tgcc.physics.analysis import analyze_llrf_nyquist, calculate_d_mode_robinson
from stable_tgcc.physics.machine_init import calc_bunch_form_factor_eq, compute_optimal_mc_detune


# Exact 7-stop pastel YlOrRd colormap from MATLAB STABLE_App.m (lines 4759-4771)
_MATLAB_YLORRD_CMAP = LinearSegmentedColormap.from_list(
    "matlab_stable_ylorrd",
    [
        (0.00, (0.98, 0.99, 0.97)),  # lowest/damped: near-white soft mint-ivory
        (0.08, (1.00, 0.98, 0.90)),  # stable: very light warm cream
        (0.18, (1.00, 0.93, 0.74)),  # boundary (~40 s^-1): soft pastel yellow
        (0.32, (1.00, 0.84, 0.56)),  # weak instability (~100 s^-1): soft apricot
        (0.50, (0.99, 0.70, 0.44)),  # moderate instability (~200 s^-1): warm salmon-orange
        (0.72, (0.96, 0.52, 0.40)),  # strong instability (~350 s^-1): soft coral red
        (1.00, (0.90, 0.34, 0.32)),  # extreme blowup (~600 s^-1): clear rose-red
    ],
    N=256,
)
_MATLAB_YLORRD_CMAP.set_bad(color=(0.95, 0.96, 0.97, 0.0))


class _UiSignals(QObject):
    telemetry_received = Signal(object)
    scan_progress_received = Signal(object)
    scan_point_received = Signal(object)
    scan_completed = Signal(int)
    job_error = Signal(str)
    log_message = Signal(str)
    status_text = Signal(str)
    async_task_step = Signal(str, int)
    connect_finished = Signal(object, str, bool, object)
    connect_failed = Signal(str, bool)
    deploy_finished = Signal(str)
    deploy_failed = Signal(str)


class STABLEMainWindow(QMainWindow):
    """Main Desktop Client Window for STABLE_Py_TGCC."""

    def __init__(
        self,
        vault: UserVaultManager,
        session: UserSession,
        startup_cb: Optional[Callable[[int, str], None]] = None,
    ):
        super().__init__()
        self.vault = vault
        self.session = session
        self.connector = STABLEClientConnector()
        self.signals = _UiSignals()
        self._busy_dialog: Optional[AsyncTaskProgressDialog] = None
        self._updating_ui_batch: bool = False

        # Caching & UI Redraw Throttling to keep GUI responsive at 60 FPS
        self._ref_curves_cache: Dict[Tuple[Any, ...], Tuple[np.ndarray, List[float], List[float], List[float], List[float]]] = {}
        self._active_marker_artist: Any = None
        self._last_waveform_draw_ts: float = 0.0
        self._last_progress_ui_ts: float = 0.0
        self._line_hc_amp: Any = None
        self._line_hc_pha: Any = None
        self._line_mc_amp: Any = None
        self._line_mc_pha: Any = None

        self._map_redraw_timer = QTimer(self)
        self._map_redraw_timer.setSingleShot(True)
        self._map_redraw_timer.setInterval(110)
        self._map_redraw_timer.timeout.connect(self._redraw_stability_map)

        self.signals.telemetry_received.connect(self._on_telemetry_frame)
        self.signals.scan_progress_received.connect(self._on_scan_progress)
        self.signals.scan_point_received.connect(self._on_scan_point)
        self.signals.scan_completed.connect(self._on_scan_completed)
        self.signals.job_error.connect(self._on_job_error)
        self.signals.log_message.connect(self._append_log)
        self.signals.status_text.connect(self._set_status)
        self.signals.async_task_step.connect(self._on_async_task_step)
        self.signals.connect_finished.connect(self._on_connect_finished)
        self.signals.connect_failed.connect(self._on_connect_failed)
        self.signals.deploy_finished.connect(self._on_deploy_finished)
        self.signals.deploy_failed.connect(self._on_deploy_failed)

        self.setWindowTitle(f"STABLE_Py_TGCC — Multi-A100 GPU RF & Beam Stability Suite [User: {session.username}]")
        self.setMinimumSize(1020, 640)
        screen = QGuiApplication.primaryScreen()
        if screen is not None:
            avail = screen.availableGeometry()
            target_w = min(1440, max(1020, int(avail.width() * 0.94)))
            target_h = min(860, max(640, int(avail.height() * 0.90)))
            self.resize(target_w, target_h)
        else:
            self.resize(1380, 820)

        # 2D Scan state matrices matching MATLAB LastScanData
        self._df_grid_kHz: np.ndarray = np.array([])
        self._i0_grid_mA: np.ndarray = np.array([])
        self._stability_code_mat: np.ndarray = np.array([])
        self._growth_rate_mat: np.ndarray = np.array([])
        self._bunch_len_mat: np.ndarray = np.array([])
        self._jitter_mat: np.ndarray = np.array([])
        self._active_scan_marker: Optional[Tuple[float, float]] = None  # (Df_kHz, I0_mA)
        self._active_scan_markers: Dict[int, Tuple[float, float]] = {}  # point_index -> (Df_kHz, I0_mA)
        self._completed_count: int = 0
        self._total_scan_points: int = 0
        self._current_job_id: str = ""
        self._map_colorbar: Any = None
        self._last_telemetry_frame: Optional[TelemetryFrame] = None
        self._offline_analysis_data: Optional[Dict[str, Any]] = None
        self._show_critical_boundary_dots: bool = False

        if startup_cb:
            startup_cb(20, "[1/4] Building parameter panels & interactive Matplotlib canvases...")
        self._build_ui()

        if startup_cb:
            startup_cb(55, "[2/4] Loading encrypted user vault & SOLEIL II presets...")
        self._load_user_profile_and_presets()

        if startup_cb:
            startup_cb(85, "[3/4] Computing initial LLRF Nyquist & Bode stability margins...")
        self._refresh_nyquist_plot()
        self._auto_load_latest_scan_on_startup()

    def post_show_auto_connect(self) -> None:
        """Called immediately after `win.show()` so backend/SSH connection never blocks startup."""
        self._connect_backend(silent=True)

    # =========================================================================
    # UI Construction
    # =========================================================================
    def _build_ui(self) -> None:
        central = QWidget()
        self.setCentralWidget(central)
        root_layout = QVBoxLayout(central)
        root_layout.setContentsMargins(8, 8, 8, 8)
        root_layout.setSpacing(6)

        # 1. Top Toolbar: User, Presets, Connection, Run/Stop Controls
        top_bar = QHBoxLayout()
        lbl_user = QLabel(f"Logged in: {self.session.username}")
        lbl_user.setFont(QFont("Segoe UI", 10, QFont.Weight.Bold))
        lbl_user.setStyleSheet("color: #1d4ed8; padding-right: 8px;")
        top_bar.addWidget(lbl_user)

        top_bar.addWidget(QLabel("Preset:"))
        self.preset_combo = QComboBox()
        self.preset_combo.setMinimumWidth(190)
        self.preset_combo.currentTextChanged.connect(self._on_preset_selected)
        top_bar.addWidget(self.preset_combo)

        btn_save_preset = QPushButton("Save Preset")
        btn_save_preset.clicked.connect(self._on_save_preset)
        top_bar.addWidget(btn_save_preset)

        btn_del_preset = QPushButton("Delete Preset")
        btn_del_preset.clicked.connect(self._on_delete_preset)
        top_bar.addWidget(btn_del_preset)

        btn_load_cfg = QPushButton("Load .mat/.m")
        btn_load_cfg.setToolTip("Load STABLE configuration (.mat) or machine script (.m)")
        btn_load_cfg.clicked.connect(self._on_load_config_file)
        top_bar.addWidget(btn_load_cfg)

        btn_save_cfg = QPushButton("Save .mat")
        btn_save_cfg.setToolTip("Save current UI configuration to a MATLAB .mat file")
        btn_save_cfg.clicked.connect(self._on_save_config_file)
        top_bar.addWidget(btn_save_cfg)

        top_bar.addSpacing(15)
        self.top_spinner = SpinningRingWidget(size=22, parent=self)
        self.top_spinner.stop()
        self.top_spinner.setVisible(False)
        top_bar.addWidget(self.top_spinner)

        self.conn_badge = QLabel("Server: Initializing...")
        self.conn_badge.setMinimumWidth(140)
        self.conn_badge.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        self.conn_badge.setStyleSheet("background: #e5e7eb; padding: 4px 8px; border-radius: 4px; font-weight: bold;")
        top_bar.addWidget(self.conn_badge, stretch=1)

        self.btn_run_single = QPushButton("▶ Run Single Tracking")
        self.btn_run_single.setStyleSheet("background-color: #16a34a; color: white; font-weight: bold; padding: 6px 14px;")
        self.btn_run_single.clicked.connect(self._start_single_tracking)
        top_bar.addWidget(self.btn_run_single)

        self.btn_run_scan = QPushButton("⊞ Run 2D Scan (Multi-A100)")
        self.btn_run_scan.setStyleSheet("background-color: #2563eb; color: white; font-weight: bold; padding: 6px 14px;")
        self.btn_run_scan.clicked.connect(self._start_2d_scan)
        top_bar.addWidget(self.btn_run_scan)

        self.btn_stop = QPushButton("■ Stop Job")
        self.btn_stop.setStyleSheet("background-color: #dc2626; color: white; font-weight: bold; padding: 6px 12px;")
        self.btn_stop.clicked.connect(self._stop_active_job)
        top_bar.addWidget(self.btn_stop)

        root_layout.addLayout(top_bar)

        # 2. Main Splitter: Left Config Tabs | Right Interactive Plot Tabs
        splitter = QSplitter(Qt.Orientation.Horizontal)
        root_layout.addWidget(splitter, stretch=1)

        self.left_tabs = QTabWidget()
        self.left_tabs.setMinimumWidth(420)
        self.left_tabs.setMaximumWidth(530)
        self.left_tabs.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Expanding)
        splitter.addWidget(self.left_tabs)

        self._build_tab_lattice()
        self._build_tab_cavities()
        self._build_tab_tracking()
        self._build_tab_feedback()
        self._build_tab_homs()
        self._build_tab_scan_and_tgcc()
        for i in range(self.left_tabs.count()):
            page = self.left_tabs.widget(i)
            page.setMinimumSize(100, 100)
            page.minimumSizeHint = lambda: QSize(380, 300)

        self.right_tabs = QTabWidget()
        self.right_tabs.setMinimumWidth(480)
        self.right_tabs.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Expanding)
        splitter.addWidget(self.right_tabs)

        self._build_plot_tab_live_tracking()
        self._build_plot_tab_stability_map()
        self._build_plot_tab_offline_analysis()
        self._build_plot_tab_nyquist()
        self._build_plot_tab_diagnostics()
        for i in range(self.right_tabs.count()):
            rpage = self.right_tabs.widget(i)
            rpage.setMinimumSize(200, 200)
            rpage.minimumSizeHint = lambda: QSize(460, 360)

        # 3. Bottom Status & Progress Bar
        bottom_bar = QHBoxLayout()
        self.status_label = QLabel("Ready.")
        self.status_label.setFont(QFont("Segoe UI", 9))
        self.status_label.setMinimumWidth(150)
        self.status_label.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        bottom_bar.addWidget(self.status_label, stretch=1)

        self.progress_bar = QProgressBar()
        self.progress_bar.setMaximumWidth(280)
        self.progress_bar.setValue(0)
        bottom_bar.addWidget(self.progress_bar)
        root_layout.addLayout(bottom_bar)

    def _make_dspin(self, val: float, lo: float, hi: float, dec: int = 4, step: float = 1.0) -> QDoubleSpinBox:
        sp = QDoubleSpinBox()
        sp.setDecimals(dec)
        sp.setRange(lo, hi)
        sp.setSingleStep(step)
        sp.setValue(val)
        return sp

    def _build_tab_lattice(self) -> None:
        w = QWidget()
        form = QFormLayout(w)
        self.f_C = self._make_dspin(354.097, 1.0, 100000.0, 3, 1.0)
        self.f_E0_GeV = self._make_dspin(2.75, 0.1, 100.0, 4, 0.05)
        self.f_I0_mA = self._make_dspin(338.0, 0.01, 5000.0, 2, 10.0)
        self.f_U0_keV = self._make_dspin(470.0, 1.0, 50000.0, 2, 10.0)
        self.f_h = QSpinBox()
        self.f_h.setRange(1, 100000)
        self.f_h.setValue(416)
        self.f_alpha_c = self._make_dspin(1.07e-4, 1e-7, 1.0, 7, 1e-5)
        self.f_tau_s_ms = self._make_dspin(13.7, 0.1, 10000.0, 3, 0.5)
        self.f_tau_z_ms = self._make_dspin(13.7, 0.1, 10000.0, 3, 0.5)
        self.f_sigma_e0 = self._make_dspin(1.06e-3, 1e-6, 0.1, 6, 1e-4)
        self.lbl_sigma_t0 = QLabel("9.00 ps")
        self.lbl_sigma_t0.setStyleSheet("font-weight: bold; color: #047857;")

        form.addRow("Circumference C [m]:", self.f_C)
        form.addRow("Beam Energy E0 [GeV]:", self.f_E0_GeV)
        form.addRow("Beam Current I0 [mA]:", self.f_I0_mA)
        form.addRow("Radiation Loss U0 [keV/turn]:", self.f_U0_keV)
        form.addRow("Harmonic Number h:", self.f_h)
        form.addRow("Momentum Compaction alpha_c:", self.f_alpha_c)
        form.addRow("Longitudinal Damping tau_s [ms]:", self.f_tau_s_ms)
        form.addRow("Energy Damping tau_z [ms]:", self.f_tau_z_ms)
        form.addRow("Natural Energy Spread sigma_e0:", self.f_sigma_e0)
        form.addRow("Natural Bunch Length sigma_t0:", self.lbl_sigma_t0)

        self.left_tabs.addTab(w, "1. Lattice")

    def _build_tab_cavities(self) -> None:
        w = QWidget()
        layout = QVBoxLayout(w)

        g_mc = QGroupBox("Main Cavity (MC)")
        f_mc = QFormLayout(g_mc)
        self.f_Ncav_mc = QSpinBox()
        self.f_Ncav_mc.setRange(1, 64)
        self.f_Ncav_mc.setValue(4)
        self.f_Q0_mc = self._make_dspin(35700.0, 10.0, 1e10, 1, 1000.0)
        self.f_RoQ_mc = self._make_dspin(140.0, 1.0, 10000.0, 2, 5.0)
        self.f_beta_mc = self._make_dspin(4.95, 0.0, 1000.0, 3, 0.1)
        self.f_V_mc_MV = self._make_dspin(1.70, 0.01, 100.0, 3, 0.05)
        self.f_Q_mc = self._make_dspin(6000.0, 10.0, 1e10, 1, 500.0)
        self.f_R_mc_MOhm = self._make_dspin(3.36, 0.001, 1e5, 4, 0.1)
        self.f_df_mc_kHz = self._make_dspin(-24.8, -5000.0, 5000.0, 3, 1.0)
        self.chk_mc_auto = QCheckBox("Robinson Optimal Detuning (Auto)")
        self.chk_mc_auto.setChecked(True)

        for sp in (self.f_Ncav_mc, self.f_Q0_mc, self.f_RoQ_mc, self.f_beta_mc):
            sp.valueChanged.connect(self._sync_cavity_derived_ui)

        f_mc.addRow("Number of Cavities Ncav_mc:", self.f_Ncav_mc)
        f_mc.addRow("Unloaded Q0_mc:", self.f_Q0_mc)
        f_mc.addRow("R/Q per Cavity [Ohm]:", self.f_RoQ_mc)
        f_mc.addRow("Coupling Factor beta_mc:", self.f_beta_mc)
        f_mc.addRow("Total Voltage V_mc [MV]:", self.f_V_mc_MV)
        f_mc.addRow("Loaded Q_mc:", self.f_Q_mc)
        f_mc.addRow("Total Loaded R_mc [MOhm]:", self.f_R_mc_MOhm)
        f_mc.addRow("MC Detuning Df_mc [kHz]:", self.f_df_mc_kHz)
        f_mc.addRow("", self.chk_mc_auto)
        layout.addWidget(g_mc)

        g_hc = QGroupBox("Passive Harmonic Cavity (HHC)")
        f_hc = QFormLayout(g_hc)
        self.f_n_hc = QSpinBox()
        self.f_n_hc.setRange(2, 16)
        self.f_n_hc.setValue(4)
        self.f_Ncav_hc = QSpinBox()
        self.f_Ncav_hc.setRange(1, 32)
        self.f_Ncav_hc.setValue(2)
        self.f_Q0_hc = self._make_dspin(31000.0, 10.0, 1e10, 1, 1000.0)
        self.f_RoQ_hc = self._make_dspin(30.0, 1.0, 10000.0, 2, 2.0)
        self.f_beta_hc = self._make_dspin(0.0, 0.0, 1000.0, 3, 0.1)
        self.f_Q_hc = self._make_dspin(31000.0, 10.0, 1e10, 1, 1000.0)
        self.f_R_hc_MOhm = self._make_dspin(1.86, 0.001, 1e6, 4, 0.1)
        self.f_df_hc_kHz = self._make_dspin(67.88, 0.1, 5000.0, 3, 2.0)
        self.chk_hc_auto = QCheckBox("Auto Flat-Potential Detuning (xi = 1.0)")
        self.chk_hc_auto.setChecked(False)
        self.chk_ff_iter = QCheckBox("Fokker-Planck Form Factor F_n(sigma_t) Iteration")
        self.chk_ff_iter.setChecked(True)
        self.chk_ff_iter.toggled.connect(lambda val: self._on_ff_iter_toggled(val, "cav"))

        btn_opt_hc = QPushButton("⚡ Optimize HC Detuning (Physical Stability Constraint, xi=0.98)")
        btn_opt_hc.setStyleSheet("background-color: #0284c7; color: white; font-weight: bold; padding: 4px 8px;")
        btn_opt_hc.clicked.connect(self._on_optimize_hc_detuning)

        for sp in (self.f_Ncav_hc, self.f_Q0_hc, self.f_RoQ_hc, self.f_beta_hc):
            sp.valueChanged.connect(self._sync_cavity_derived_ui)

        f_hc.addRow("Harmonic Multiple n_hc:", self.f_n_hc)
        f_hc.addRow("Number of Cavities Ncav_hc:", self.f_Ncav_hc)
        f_hc.addRow("Unloaded Q0_hc:", self.f_Q0_hc)
        f_hc.addRow("R/Q per Cavity [Ohm]:", self.f_RoQ_hc)
        f_hc.addRow("Coupling beta_hc:", self.f_beta_hc)
        f_hc.addRow("Loaded Q_hc:", self.f_Q_hc)
        f_hc.addRow("Total Loaded R_hc [MOhm]:", self.f_R_hc_MOhm)
        f_hc.addRow("HHC Detuning Df_hc [kHz]:", self.f_df_hc_kHz)
        f_hc.addRow("", self.chk_hc_auto)
        f_hc.addRow("", self.chk_ff_iter)
        f_hc.addRow("", btn_opt_hc)
        layout.addWidget(g_hc)

        self.left_tabs.addTab(w, "2. Cavities")

    def _build_tab_tracking(self) -> None:
        w = QWidget()
        form = QFormLayout(w)
        self.f_pattern = QComboBox()
        self.f_pattern.addItems(["Uniform (100%)", "80% filling", "90% filling", "Single bunch", "Custom"])
        self.f_custom_pat = QLineEdit()
        self.f_custom_pat.setPlaceholderText("e.g., np.ones(h)")
        self.f_Par_num = QSpinBox()
        self.f_Par_num.setRange(50, 1000000)
        self.f_Par_num.setSingleStep(1000)
        self.f_Par_num.setValue(10000)
        self.f_Track_num = QSpinBox()
        self.f_Track_num.setRange(50, 2000000)
        self.f_Track_num.setSingleStep(5000)
        self.f_Track_num.setValue(30000)
        self.f_Recor_step = QSpinBox()
        self.f_Recor_step.setRange(1, 5000)
        self.f_Recor_step.setValue(10)
        self.chk_short_wake = QCheckBox("Enable Intra-Bunch Short-Range Wake")
        self.chk_short_wake.setChecked(True)

        form.addRow("Filling Pattern:", self.f_pattern)
        form.addRow("Custom Pattern Expr:", self.f_custom_pat)
        form.addRow("Macro-Particles / Bunch:", self.f_Par_num)
        form.addRow("Tracking Turns:", self.f_Track_num)
        form.addRow("Record Interval [turns]:", self.f_Recor_step)
        form.addRow("", self.chk_short_wake)

        self.left_tabs.addTab(w, "3. Tracking")

    def _build_tab_feedback(self) -> None:
        w = QWidget()
        layout = QVBoxLayout(w)

        g_fb_pre = QGroupBox("Feedback Quick Presets (ξ=1 Stability Island & LLRF Tuning)")
        v_pre = QVBoxLayout(g_fb_pre)
        self.fb_preset_combo = QComboBox()
        self.fb_preset_combo.addItems([
            "SOLEIL II Default (DFB Auto 0.25-0.20, PI kp=0.8, ki=1e-4)",
            "mbtrack2 Weak DFB (G_dfb=0.15, PI kp=0.8, ki=1e-4)",
            "mbtrack2 Moderate DFB (G_dfb=0.25, PI kp=1.0, ki=1e-4)",
            "PI Only (DFB Off, kp=0.8, ki=1e-4)",
            "Open Loop (All Feedback Off)",
        ])
        v_pre.addWidget(self.fb_preset_combo)
        h_pre_btn = QHBoxLayout()
        btn_apply_fb = QPushButton("Apply Feedback Preset")
        btn_apply_fb.setStyleSheet("color: #0369a1; font-weight: bold;")
        btn_apply_fb.clicked.connect(self._on_apply_feedback_preset)
        btn_save_soleil = QPushButton("Save to SOLEIL II Default")
        btn_save_soleil.setStyleSheet("background-color: #059669; color: white; font-weight: bold; padding: 3px 8px;")
        btn_save_soleil.setToolTip("Save all current UI settings into the 'SOLEIL II (Default)' preset")
        btn_save_soleil.clicked.connect(self._on_save_to_soleil_default)
        h_pre_btn.addWidget(btn_apply_fb)
        h_pre_btn.addWidget(btn_save_soleil)
        v_pre.addLayout(h_pre_btn)
        layout.addWidget(g_fb_pre)

        g_pi = QGroupBox("Proportional-Integral (PI) Cavity Loop")
        f_pi = QFormLayout(g_pi)
        self.chk_pi_enable = QCheckBox("Enable PI Feedback")
        self.chk_pi_enable.setChecked(True)
        self.f_pi_kp = self._make_dspin(0.8, 0.0, 50.0, 4, 0.1)
        self.f_pi_ki = self._make_dspin(1.0e-4, 0.0, 1.0, 7, 1e-5)
        self.f_pi_m = QSpinBox()
        self.f_pi_m.setRange(1, 10000)
        self.f_pi_m.setValue(80)
        self.f_pi_d = QSpinBox()
        self.f_pi_d.setRange(1, 50000)
        self.f_pi_d.setValue(800)

        f_pi.addRow("", self.chk_pi_enable)
        f_pi.addRow("Proportional Gain KP:", self.f_pi_kp)
        f_pi.addRow("Integral Gain KI:", self.f_pi_ki)
        f_pi.addRow("Sample Window PI.m [buckets]:", self.f_pi_m)
        f_pi.addRow("Loop Delay PI.d [buckets]:", self.f_pi_d)
        layout.addWidget(g_pi)

        g_dfb = QGroupBox("Direct RF Feedback (DFB) Loop (Phase Shift = 0° is Neg FB)")
        f_dfb = QFormLayout(g_dfb)
        self.chk_dfb_enable = QCheckBox("Enable Direct RF Feedback (DFB)")
        self.chk_dfb_enable.setChecked(True)
        self.f_dfb_gain = self._make_dspin(0.20, 0.0, 5.0, 3, 0.05)
        self.f_dfb_phase = self._make_dspin(0.0, -180.0, 180.0, 2, 1.0)
        self.f_dfb_delay = QSpinBox()
        self.f_dfb_delay.setRange(1, 50000)
        self.f_dfb_delay.setValue(704)
        self.f_dfb_m = QSpinBox()
        self.f_dfb_m.setRange(1, 10000)
        self.f_dfb_m.setValue(208)
        self.f_dfb_every = QSpinBox()
        self.f_dfb_every.setRange(1, 10000)
        self.f_dfb_every.setValue(208)
        self.chk_dfb_autoscale = QCheckBox("Auto-Scale Gain with Current (I1..I2 -> G1..G2)")
        self.chk_dfb_autoscale.setChecked(True)
        self.chk_dfb_adaptive = QCheckBox("Adaptive Margin-Constrained Gain")
        self.chk_dfb_adaptive.setChecked(False)
        self.f_dfb_I1 = self._make_dspin(300.0, 1.0, 5000.0, 1, 25.0)
        self.f_dfb_G1 = self._make_dspin(0.25, 0.0, 5.0, 3, 0.02)
        self.f_dfb_I2 = self._make_dspin(500.0, 1.0, 5000.0, 1, 25.0)
        self.f_dfb_G2 = self._make_dspin(0.20, 0.0, 5.0, 3, 0.02)

        btn_upd_nyq = QPushButton("Update LLRF Nyquist & Bode")
        btn_upd_nyq.clicked.connect(self._refresh_nyquist_plot)

        f_dfb.addRow("", self.chk_dfb_enable)
        f_dfb.addRow("DFB Gain:", self.f_dfb_gain)
        f_dfb.addRow("Phase Shift [deg] (0° = Neg FB):", self.f_dfb_phase)
        f_dfb.addRow("Loop Delay [buckets]:", self.f_dfb_delay)
        f_dfb.addRow("Sample Window m [buckets]:", self.f_dfb_m)
        f_dfb.addRow("Update Interval every [buckets]:", self.f_dfb_every)
        f_dfb.addRow("", self.chk_dfb_autoscale)
        f_dfb.addRow("", self.chk_dfb_adaptive)
        f_dfb.addRow("Anchor I1 [mA] / Gain G1:", self._pair_widget(self.f_dfb_I1, self.f_dfb_G1))
        f_dfb.addRow("Anchor I2 [mA] / Gain G2:", self._pair_widget(self.f_dfb_I2, self.f_dfb_G2))
        f_dfb.addRow("", btn_upd_nyq)
        layout.addWidget(g_dfb)

        self.left_tabs.addTab(w, "4. Feedback")

    def _pair_widget(self, w1: QWidget, w2: QWidget) -> QWidget:
        p = QWidget()
        h = QHBoxLayout(p)
        h.setContentsMargins(0, 0, 0, 0)
        h.addWidget(w1)
        h.addWidget(w2)
        return p

    def _build_tab_homs(self) -> None:
        w = QWidget()
        layout = QVBoxLayout(w)

        self.chk_hom_align = QCheckBox("Worst-Case Synchrotron Sideband Alignment (Auto-Align)")
        self.chk_hom_align.setChecked(True)
        layout.addWidget(self.chk_hom_align)

        self.hom_table = QTableWidget(0, 3)
        self.hom_table.setHorizontalHeaderLabels(["Q_hom", "R_hom [Ohm]", "Freq [Hz]"])
        self.hom_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        layout.addWidget(self.hom_table, stretch=1)

        btn_row = QHBoxLayout()
        btn_add = QPushButton("+ Add HOM Mode")
        btn_add.clicked.connect(lambda: self._append_hom_row(830.0, 6500.0, 1.7e9))
        btn_del = QPushButton("- Remove Selected")
        btn_del.clicked.connect(self._remove_selected_hom)
        btn_align = QPushButton("Align Sideband")
        btn_align.setStyleSheet("color: #0369a1; font-weight: bold;")
        btn_align.setToolTip("Align HOM frequencies to worst-case upper synchrotron sideband (mbtrack2 standard)")
        btn_align.clicked.connect(self._on_align_homs_button)
        btn_clear = QPushButton("Clear All")
        btn_clear.clicked.connect(lambda: self.hom_table.setRowCount(0))
        btn_row.addWidget(btn_add)
        btn_row.addWidget(btn_del)
        btn_row.addWidget(btn_align)
        btn_row.addWidget(btn_clear)
        layout.addLayout(btn_row)

        self.left_tabs.addTab(w, "5. HOMs")

    def _append_hom_row(self, q_val: float, r_val: float, f_val: float) -> None:
        r = self.hom_table.rowCount()
        self.hom_table.insertRow(r)
        self.hom_table.setItem(r, 0, QTableWidgetItem(f"{q_val:.1f}"))
        self.hom_table.setItem(r, 1, QTableWidgetItem(f"{r_val:.1f}"))
        self.hom_table.setItem(r, 2, QTableWidgetItem(f"{f_val:.6e}"))

    def _remove_selected_hom(self) -> None:
        row = self.hom_table.currentRow()
        if row >= 0:
            self.hom_table.removeRow(row)

    def _build_tab_scan_and_tgcc(self) -> None:
        w = QWidget()
        layout = QVBoxLayout(w)

        g_scan = QGroupBox("2D Stability Map Scan Grid & Multi-A100 Dispatch")
        f_scan = QFormLayout(g_scan)
        self.f_scan_df_min = self._make_dspin(50.0, 1.0, 5000.0, 1, 5.0)
        self.f_scan_df_max = self._make_dspin(150.0, 1.0, 5000.0, 1, 5.0)
        self.f_scan_df_step = self._make_dspin(5.0, 0.5, 500.0, 1, 1.0)
        self.f_scan_i0_min = self._make_dspin(300.0, 1.0, 5000.0, 1, 25.0)
        self.f_scan_i0_max = self._make_dspin(500.0, 1.0, 5000.0, 1, 25.0)
        self.f_scan_i0_step = self._make_dspin(25.0, 1.0, 1000.0, 1, 5.0)
        self.chk_dmode_prescreen = QCheckBox("D-Mode Robinson Pre-Screening (Superconducting Q_hc >= 1e5 only)")
        self.chk_dmode_prescreen.setChecked(False)
        self.chk_dmode_prescreen.setToolTip(
            "Analytical D-mode Robinson pre-filter (T. He PRAB 2024). Only valid for superconducting high-Q cavities (Q_hc >= 1e5); normal-conducting SOLEIL II cavities always use full particle tracking."
        )
        self.f_Q_hc.valueChanged.connect(lambda _v: self._sync_dmode_prescreen_availability())
        self.f_scan_gpus = QSpinBox()
        self.f_scan_gpus.setRange(1, 16)
        self.f_scan_gpus.setValue(1)
        self.f_scan_gpus.setToolTip(
            "Number of A100 GPUs to request on Topaze (1 or 2 GPUs start immediately without queue wait)."
        )
        self.f_workers_per_gpu = QSpinBox()
        self.f_workers_per_gpu.setRange(1, 32)
        self.f_workers_per_gpu.setValue(16)
        self.f_workers_per_gpu.setToolTip(
            "Concurrent grid points per A100-80GB GPU. Sweet spot: 12~16 pts/GPU (~5.6 GB / 80 GB VRAM, 16 of 32 CPU cores). Max: 32."
        )
        self.lbl_gpu_concurrency = QLabel("1 A100 × 16 Pts/GPU = 16 points in parallel (~5.6 GB / 80 GB VRAM)")
        self.lbl_gpu_concurrency.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        self.lbl_gpu_concurrency.setStyleSheet("color: #0369a1; font-weight: bold;")
        self.f_exec_mode = QComboBox()
        self.f_exec_mode.addItems(["ccc_msub_batch", "daemon_stream"])
        self.f_walltime_hours = self._make_dspin(2.5, 0.5, 24.0, 1, 0.5)
        self.chk_auto_walltime = QCheckBox("Auto-Estimate")
        self.chk_auto_walltime.setChecked(True)
        self.chk_auto_walltime.toggled.connect(self._update_estimated_walltime_ui)
        for w_sp in (
            self.f_scan_df_min,
            self.f_scan_df_max,
            self.f_scan_df_step,
            self.f_scan_i0_min,
            self.f_scan_i0_max,
            self.f_scan_i0_step,
            self.f_scan_gpus,
            self.f_workers_per_gpu,
            self.f_Par_num,
            self.f_Track_num,
        ):
            w_sp.valueChanged.connect(self._update_estimated_walltime_ui)

        f_scan.addRow("Df_HC Min / Max [kHz]:", self._pair_widget(self.f_scan_df_min, self.f_scan_df_max))
        f_scan.addRow("Df_HC Step [kHz]:", self.f_scan_df_step)
        f_scan.addRow("I0 Min / Max [mA]:", self._pair_widget(self.f_scan_i0_min, self.f_scan_i0_max))
        f_scan.addRow("I0 Step [mA]:", self.f_scan_i0_step)
        f_scan.addRow("Requested A100 GPUs:", self.f_scan_gpus)
        f_scan.addRow("Points / A100 GPU (单卡并发):", self.f_workers_per_gpu)
        f_scan.addRow("Total Concurrency:", self.lbl_gpu_concurrency)
        f_scan.addRow("Job Walltime [hours]:", self._pair_widget(self.f_walltime_hours, self.chk_auto_walltime))
        f_scan.addRow("TGCC Execution Mode:", self.f_exec_mode)
        f_scan.addRow("", self.chk_dmode_prescreen)
        layout.addWidget(g_scan)

        g_tgcc = QGroupBox("TGCC Topaze A100 SSH Connection (AES-Encrypted in User Vault)")
        f_tgcc = QFormLayout(g_tgcc)
        self.f_conn_mode = QComboBox()
        self.f_conn_mode.addItems(["local", "ssh_tgcc"])
        self.f_ssh_host = QLineEdit("topaze.ccc.cea.fr")
        self.f_ssh_port = QSpinBox()
        self.f_ssh_port.setRange(1, 65535)
        self.f_ssh_port.setValue(22)
        self.f_ssh_user = QLineEdit(self.session.username or "zhaolu")
        self.f_ssh_pwd = QLineEdit()
        self.f_ssh_pwd.setEchoMode(QLineEdit.EchoMode.Password)
        self.f_ssh_pwd.setPlaceholderText("Encrypted SSH Password / Passphrase")
        self.f_ssh_otp = QLineEdit()
        self.f_ssh_otp.setPlaceholderText("Optional 2FA / OTP token for TGCC login")
        self.f_ccc_proj = QLineEdit("soldai")
        self.f_ccc_part = QLineEdit("a100")
        self.f_ccc_qos = QComboBox()
        self.f_ccc_qos.addItems(["normal", "test", "long"])
        self.f_ccc_qos.setToolTip(
            "TGCC QoS (#MSUB -Q, see ccc_mqinfo): 'test' gives high-priority boost for short jobs (<=30 min); 'normal' for standard runs."
        )
        self.f_remote_dir = QLineEdit("/ccc/work/cont003/soleil/zhaolu/STABLE_Py_TGCC")
        self.f_python_cmd = QLineEdit("python3")

        f_tgcc.addRow("Server Mode:", self.f_conn_mode)
        f_tgcc.addRow("TGCC SSH Host / Port:", self._pair_widget(self.f_ssh_host, self.f_ssh_port))
        f_tgcc.addRow("SSH Username:", self.f_ssh_user)
        f_tgcc.addRow("SSH Password (Encrypted):", self.f_ssh_pwd)
        f_tgcc.addRow("2FA / OTP Code (Session):", self.f_ssh_otp)
        f_tgcc.addRow("TGCC Account (-A) / Queue (-q):", self._pair_widget(self.f_ccc_proj, self.f_ccc_part))
        f_tgcc.addRow("TGCC Priority QoS (-Q):", self.f_ccc_qos)
        f_tgcc.addRow("Remote Work Dir:", self.f_remote_dir)
        f_tgcc.addRow("Remote Python Cmd:", self.f_python_cmd)

        btn_row = QHBoxLayout()
        btn_save_conn = QPushButton("Save SSH Profile")
        btn_save_conn.clicked.connect(self._on_save_tgcc_profile)
        btn_connect = QPushButton("Connect / Test Server")
        btn_connect.setStyleSheet("background-color: #0284c7; color: white; font-weight: bold;")
        btn_connect.clicked.connect(lambda: self._connect_backend(silent=False))
        btn_deploy = QPushButton("Deploy Code")
        btn_deploy.clicked.connect(self._on_deploy_to_tgcc)
        btn_clean_q = QPushButton("Clean Pending Queue")
        btn_clean_q.setStyleSheet("color: #b91c1c; font-weight: bold;")
        btn_clean_q.setToolTip("Cancel any stuck/pending STABLE_* jobs in Topaze queue and show SLURM reason.")
        btn_clean_q.clicked.connect(self._on_clean_pending_queue)
        btn_row.addWidget(btn_save_conn)
        btn_row.addWidget(btn_connect)
        btn_row.addWidget(btn_deploy)
        btn_row.addWidget(btn_clean_q)
        f_tgcc.addRow("", btn_row)
        layout.addWidget(g_tgcc)

        self.left_tabs.addTab(w, "6. Scan & TGCC")

    # =========================================================================
    # Right Interactive Plot Tabs (100% Matched to MATLAB STABLE_App.m)
    # =========================================================================
    def _build_plot_tab_live_tracking(self) -> None:
        w = QWidget()
        v = QVBoxLayout(w)
        top_ctrl = QHBoxLayout()
        top_ctrl.addWidget(QLabel("Tracking View:"))
        self.live_view_combo = QComboBox()
        self.live_view_combo.addItems([
            "Cavity Voltages |V| & Phase vs Bucket (HC & MC)",
            "Centroid <tau> & Bunch Length sigma_tau vs Turns",
        ])
        self.live_view_combo.currentIndexChanged.connect(self._on_live_view_mode_changed)
        top_ctrl.addWidget(self.live_view_combo)
        top_ctrl.addStretch()

        self.fig_live = Figure(figsize=(8.0, 6.0), tight_layout=True)
        self.canvas_live = FigureCanvas(self.fig_live)
        self.canvas_live.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Ignored)
        self.canvas_live.minimumSizeHint = lambda: QSize(360, 240)
        self.nav_live = NavigationToolbar(self.canvas_live, w)
        top_ctrl.addWidget(self.nav_live)
        v.addLayout(top_ctrl)

        self.ax_hc_amp = self.fig_live.add_subplot(2, 2, 1)
        self.ax_mc_amp = self.fig_live.add_subplot(2, 2, 2)
        self.ax_hc_pha = self.fig_live.add_subplot(2, 2, 3)
        self.ax_mc_pha = self.fig_live.add_subplot(2, 2, 4)
        v.addWidget(self.canvas_live, stretch=1)
        self.right_tabs.addTab(w, "Live Tracking Waveforms")

    def _build_plot_tab_stability_map(self) -> None:
        w = QWidget()
        v = QVBoxLayout(w)
        ctrl = QHBoxLayout()
        ctrl.addWidget(QLabel("Coordinate:"))
        self.map_xaxis_combo = QComboBox()
        self.map_xaxis_combo.addItems(["HHC Detuning Df_HC [kHz]", "HHC Tuning Angle psi_HC [deg]"])
        self.map_xaxis_combo.currentIndexChanged.connect(self._redraw_stability_map)
        ctrl.addWidget(self.map_xaxis_combo)

        btn_load_mat = QPushButton("Load Scan Data...")
        btn_load_mat.setToolTip("Load scan_summary.mat and automatically restore its machine & feedback settings to UI")
        btn_load_mat.clicked.connect(self._on_load_matlab_scan_summary)
        ctrl.addWidget(btn_load_mat)

        btn_load_dir = QPushButton("Load Folder...")
        btn_load_dir.setToolTip("Select a Scan_* folder containing scan_summary.mat and SOLEIL_I0*.mat files")
        btn_load_dir.clicked.connect(self._on_load_scan_folder)
        ctrl.addWidget(btn_load_dir)

        btn_save_mat = QPushButton("Save Scan .mat...")
        btn_save_mat.setToolTip("Export current 2D scan map and settings to a MATLAB-compatible scan_summary.mat")
        btn_save_mat.clicked.connect(self._on_save_scan_summary_mat)
        ctrl.addWidget(btn_save_mat)

        btn_refresh = QPushButton("Refresh Map")
        btn_refresh.clicked.connect(self._redraw_stability_map)
        ctrl.addWidget(btn_refresh)

        btn_fit_bnd = QPushButton("Fit Safety Boundary")
        btn_fit_bnd.setStyleSheet("background-color: #d95919; color: white; font-weight: bold; padding: 4px 8px;")
        btn_fit_bnd.setToolTip("Extract critical stability boundary and fit engineering equation Df_crit(I0)")
        btn_fit_bnd.clicked.connect(self._on_fit_safety_boundary)
        ctrl.addWidget(btn_fit_bnd)

        self.chk_map_ff_iter = QCheckBox("Form Factor F_n(\u03c3_\u03c4)")
        self.chk_map_ff_iter.setChecked(True)
        self.chk_map_ff_iter.setToolTip(
            "Include self-consistent stretched bunch form factor F_n(sigma_tau) in Flat-pot. (xi=1.0) & xi=1.05 curves"
        )
        self.chk_map_ff_iter.toggled.connect(lambda val: self._on_ff_iter_toggled(val, "map"))
        ctrl.addWidget(self.chk_map_ff_iter)

        ctrl.addStretch()
        v.addLayout(ctrl)

        self.fig_map = Figure(figsize=(8.2, 5.8), tight_layout=True)
        self.canvas_map = FigureCanvas(self.fig_map)
        self.canvas_map.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Ignored)
        self.canvas_map.minimumSizeHint = lambda: QSize(360, 240)
        nav_row = QHBoxLayout()
        self.nav_map = NavigationToolbar(self.canvas_map, w)
        nav_row.addWidget(self.nav_map, stretch=1)
        nav_row.addWidget(QLabel("Legend:"))
        self.map_legend_combo = QComboBox()
        self.map_legend_combo.addItems(["Lower Right", "Upper Right", "Upper Left", "Lower Left", "Best (Auto)", "Hide"])
        self.map_legend_combo.setToolTip("Choose legend corner or drag the legend directly on the plot with the mouse")
        self.map_legend_combo.currentIndexChanged.connect(self._redraw_stability_map)
        nav_row.addWidget(self.map_legend_combo)
        v.addLayout(nav_row)
        self.ax_map = self.fig_map.add_subplot(1, 1, 1)
        v.addWidget(self.canvas_map, stretch=1)
        self.right_tabs.addTab(w, "2D Stability Map")

    def _build_plot_tab_offline_analysis(self) -> None:
        w = QWidget()
        v = QVBoxLayout(w)
        ctrl = QHBoxLayout()
        ctrl.addWidget(QLabel("Plot Mode:"))
        self.analysis_mode_combo = QComboBox()
        self.analysis_mode_combo.addItems([
            "All",
            "Bunch Length",
            "Centroid",
            "Energy Spread",
            "Bunch Length + Centroid",
            "Bunch Length + Energy Spread",
            "Centroid + Energy Spread",
            "Mode-0 Spectrum (D-mode Check)",
        ])
        self.analysis_mode_combo.currentIndexChanged.connect(self._redraw_offline_analysis_plot)
        ctrl.addWidget(self.analysis_mode_combo)

        btn_an_dir = QPushButton("Analyze Scan Folder (SOLEIL_I0*.mat)...")
        btn_an_dir.setStyleSheet("background-color: #0284c7; color: white; font-weight: bold; padding: 4px 10px;")
        btn_an_dir.clicked.connect(self._on_analyze_scan_folder)
        ctrl.addWidget(btn_an_dir)

        btn_an_mem = QPushButton("Plot Current 2D Scan Slice")
        btn_an_mem.clicked.connect(self._on_plot_current_scan_slice)
        ctrl.addWidget(btn_an_mem)
        ctrl.addStretch()

        self.fig_an = Figure(figsize=(8.0, 5.6), tight_layout=True)
        self.canvas_an = FigureCanvas(self.fig_an)
        self.canvas_an.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Ignored)
        self.canvas_an.minimumSizeHint = lambda: QSize(360, 240)
        self.nav_an = NavigationToolbar(self.canvas_an, w)
        ctrl.addWidget(self.nav_an)
        v.addLayout(ctrl)

        self.ax_an = self.fig_an.add_subplot(1, 1, 1)
        v.addWidget(self.canvas_an, stretch=1)
        self.right_tabs.addTab(w, "Analysis (BL / Centroid)")

    def _build_plot_tab_nyquist(self) -> None:
        w = QWidget()
        v = QVBoxLayout(w)
        ctrl = QHBoxLayout()
        ctrl.addWidget(QLabel("Loop Mode:"))
        self.llrf_mode_combo = QComboBox()
        self.llrf_mode_combo.addItems([
            "Gain Sweep (Compare G_dfb)",
            "Combined (DFB + PI)",
            "Direct RF Feedback (DFB Only)",
            "PI Feedback (PI Only)",
        ])
        self.llrf_mode_combo.currentIndexChanged.connect(self._refresh_nyquist_plot)
        ctrl.addWidget(self.llrf_mode_combo)

        self.chk_llrf_detune = QCheckBox("Include MC Detuning \u03c8")
        self.chk_llrf_detune.setChecked(True)
        self.chk_llrf_detune.toggled.connect(self._refresh_nyquist_plot)
        ctrl.addWidget(self.chk_llrf_detune)

        btn_upd_llrf = QPushButton("Update Plot")
        btn_upd_llrf.setStyleSheet("background-color: #2563eb; color: white; font-weight: bold; padding: 4px 10px;")
        btn_upd_llrf.clicked.connect(self._refresh_nyquist_plot)
        ctrl.addWidget(btn_upd_llrf)

        self.lbl_llrf_metrics = QLabel("GM: -- dB | PM: -- deg | Max G_dfb: --")
        self.lbl_llrf_metrics.setStyleSheet("color: #1d4ed8; font-weight: bold;")
        self.lbl_llrf_metrics.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        ctrl.addWidget(self.lbl_llrf_metrics, stretch=1)
        v.addLayout(ctrl)

        self.fig_nyq = Figure(figsize=(8.0, 5.6), tight_layout=True)
        self.canvas_nyq = FigureCanvas(self.fig_nyq)
        self.canvas_nyq.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Ignored)
        self.canvas_nyq.minimumSizeHint = lambda: QSize(360, 240)
        self.nav_nyq = NavigationToolbar(self.canvas_nyq, w)
        v.addWidget(self.nav_nyq)
        self.ax_bode_mag = self.fig_nyq.add_subplot(2, 2, 1)
        self.ax_bode_pha = self.fig_nyq.add_subplot(2, 2, 3)
        self.ax_nyquist = self.fig_nyq.add_subplot(1, 2, 2)
        v.addWidget(self.canvas_nyq, stretch=1)
        self.right_tabs.addTab(w, "LLRF Nyquist & Bode")

    def _build_plot_tab_diagnostics(self) -> None:
        w = QWidget()
        v = QVBoxLayout(w)
        self.fig_diag = Figure(figsize=(8.0, 3.8), tight_layout=True)
        self.canvas_diag = FigureCanvas(self.fig_diag)
        self.canvas_diag.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Ignored)
        self.canvas_diag.minimumSizeHint = lambda: QSize(360, 200)
        self.ax_evo = self.fig_diag.add_subplot(1, 2, 1)
        self.ax_modes = self.fig_diag.add_subplot(1, 2, 2)
        v.addWidget(self.canvas_diag, stretch=2)

        self.log_console = QTextEdit()
        self.log_console.setReadOnly(True)
        self.log_console.setFont(QFont("Consolas", 9))
        self.log_console.document().setMaximumBlockCount(1500)
        v.addWidget(self.log_console, stretch=1)
        self.right_tabs.addTab(w, "Mode Spectrum & TGCC Console")

    # =========================================================================
    # State Sync: UI <-> Dataclasses <-> User Vault
    # =========================================================================
    def _sync_cavity_derived_ui(self) -> None:
        if self._updating_ui_batch:
            return
        q_mc = self.f_Q0_mc.value() / (1.0 + max(0.0, self.f_beta_mc.value()))
        r_mc = self.f_Ncav_mc.value() * self.f_RoQ_mc.value() * q_mc
        self.f_Q_mc.setValue(q_mc)
        self.f_R_mc_MOhm.setValue(r_mc * 1e-6)

        q_hc = self.f_Q0_hc.value() / (1.0 + max(0.0, self.f_beta_hc.value()))
        r_hc = self.f_Ncav_hc.value() * self.f_RoQ_hc.value() * q_hc
        self.f_Q_hc.setValue(q_hc)
        self.f_R_hc_MOhm.setValue(r_hc * 1e-6)

    def _collect_machine_config(self) -> MachineConfig:
        homs: List[HOMMode] = []
        for r in range(self.hom_table.rowCount()):
            try:
                q_v = float(self.hom_table.item(r, 0).text())
                r_v = float(self.hom_table.item(r, 1).text())
                f_v = float(self.hom_table.item(r, 2).text())
                homs.append(HOMMode(Q=q_v, R=r_v, freq_Hz=f_v))
            except Exception:
                pass

        return MachineConfig(
            preset_name=self.preset_combo.currentText() or "Custom",
            C=self.f_C.value(),
            E0=self.f_E0_GeV.value() * 1e9,
            I0=self.f_I0_mA.value() * 1e-3,
            U0=self.f_U0_keV.value() * 1e3,
            h=self.f_h.value(),
            alpha_c=self.f_alpha_c.value(),
            tau_s=self.f_tau_s_ms.value() * 1e-3,
            tau_z=self.f_tau_z_ms.value() * 1e-3,
            sigma_e0=self.f_sigma_e0.value(),
            Ncav_mc=self.f_Ncav_mc.value(),
            Q0_mc=self.f_Q0_mc.value(),
            RoQ_mc=self.f_RoQ_mc.value(),
            beta_mc=self.f_beta_mc.value(),
            V_mc=self.f_V_mc_MV.value() * 1e6,
            Q_mc=self.f_Q_mc.value(),
            R_mc=self.f_R_mc_MOhm.value() * 1e6,
            fre_shift_mc=self.f_df_mc_kHz.value() * 1e3,
            mc_auto_detune=self.chk_mc_auto.isChecked(),
            n_hc=self.f_n_hc.value(),
            Ncav_hc=self.f_Ncav_hc.value(),
            Q0_hc=self.f_Q0_hc.value(),
            RoQ_hc=self.f_RoQ_hc.value(),
            beta_hc=self.f_beta_hc.value(),
            Q_hc=self.f_Q_hc.value(),
            R_hc=self.f_R_hc_MOhm.value() * 1e6,
            fre_shift_hc=self.f_df_hc_kHz.value() * 1e3,
            hc_auto_detune=self.chk_hc_auto.isChecked(),
            use_form_factor_iter=self.chk_ff_iter.isChecked(),
            pattern_mode=self.f_pattern.currentText(),
            custom_pattern_expr=self.f_custom_pat.text(),
            Par_num=self.f_Par_num.value(),
            Track_num=self.f_Track_num.value(),
            Recor_step=self.f_Recor_step.value(),
            short_range_wake=self.chk_short_wake.isChecked(),
            pi_enable=self.chk_pi_enable.isChecked(),
            pi_kp=self.f_pi_kp.value(),
            pi_ki=self.f_pi_ki.value(),
            pi_m=self.f_pi_m.value(),
            pi_d=self.f_pi_d.value(),
            dfb_enable=self.chk_dfb_enable.isChecked(),
            dfb_gain=self.f_dfb_gain.value(),
            dfb_phase_shift_deg=self.f_dfb_phase.value(),
            dfb_delay=self.f_dfb_delay.value(),
            dfb_m=self.f_dfb_m.value(),
            dfb_every=self.f_dfb_every.value(),
            dfb_auto_scale=self.chk_dfb_autoscale.isChecked(),
            dfb_adaptive=self.chk_dfb_adaptive.isChecked(),
            dfb_I1_mA=self.f_dfb_I1.value(),
            dfb_G1=self.f_dfb_G1.value(),
            dfb_I2_mA=self.f_dfb_I2.value(),
            dfb_G2=self.f_dfb_G2.value(),
            homs=homs,
            hom_auto_align=self.chk_hom_align.isChecked(),
        )

    def _apply_machine_config_to_ui(self, cfg: MachineConfig) -> None:
        self._updating_ui_batch = True
        try:
            self.f_C.setValue(cfg.C)
            self.f_E0_GeV.setValue(cfg.E0 * 1e-9)
            self.f_I0_mA.setValue(cfg.I0 * 1e3)
            self.f_U0_keV.setValue(cfg.U0 * 1e-3)
            self.f_h.setValue(cfg.h)
            self.f_alpha_c.setValue(cfg.alpha_c)
            self.f_tau_s_ms.setValue(cfg.tau_s * 1e3)
            self.f_tau_z_ms.setValue(cfg.tau_z * 1e3)
            self.f_sigma_e0.setValue(cfg.sigma_e0)

            self.f_Ncav_mc.setValue(cfg.Ncav_mc)
            self.f_Q0_mc.setValue(cfg.Q0_mc)
            self.f_RoQ_mc.setValue(cfg.RoQ_mc)
            self.f_beta_mc.setValue(cfg.beta_mc)
            self.f_V_mc_MV.setValue(cfg.V_mc * 1e-6)
            self.f_Q_mc.setValue(cfg.Q_mc)
            self.f_R_mc_MOhm.setValue(cfg.R_mc * 1e-6)
            self.f_df_mc_kHz.setValue(cfg.fre_shift_mc * 1e-3)
            self.chk_mc_auto.setChecked(cfg.mc_auto_detune)

            self.f_n_hc.setValue(cfg.n_hc)
            self.f_Ncav_hc.setValue(cfg.Ncav_hc)
            self.f_Q0_hc.setValue(cfg.Q0_hc)
            self.f_RoQ_hc.setValue(cfg.RoQ_hc)
            self.f_beta_hc.setValue(cfg.beta_hc)
            self.f_Q_hc.setValue(cfg.Q_hc)
            self.f_R_hc_MOhm.setValue(cfg.R_hc * 1e-6)
            self.f_df_hc_kHz.setValue(cfg.fre_shift_hc * 1e-3)
            self.chk_hc_auto.setChecked(cfg.hc_auto_detune)
            self.chk_ff_iter.setChecked(cfg.use_form_factor_iter)

            idx_p = self.f_pattern.findText(cfg.pattern_mode)
            if idx_p >= 0:
                self.f_pattern.setCurrentIndex(idx_p)
            self.f_custom_pat.setText(cfg.custom_pattern_expr)
            self.f_Par_num.setValue(cfg.Par_num)
            self.f_Track_num.setValue(cfg.Track_num)
            self.f_Recor_step.setValue(cfg.Recor_step)
            self.chk_short_wake.setChecked(cfg.short_range_wake)

            self.chk_pi_enable.setChecked(cfg.pi_enable)
            self.f_pi_kp.setValue(cfg.pi_kp)
            self.f_pi_ki.setValue(cfg.pi_ki)
            self.f_pi_m.setValue(cfg.pi_m)
            self.f_pi_d.setValue(cfg.pi_d)

            self.chk_dfb_enable.setChecked(cfg.dfb_enable)
            self.f_dfb_gain.setValue(cfg.dfb_gain)
            self.f_dfb_phase.setValue(cfg.dfb_phase_shift_deg)
            self.f_dfb_delay.setValue(cfg.dfb_delay)
            self.f_dfb_m.setValue(cfg.dfb_m)
            self.f_dfb_every.setValue(cfg.dfb_every)
            self.chk_dfb_autoscale.setChecked(cfg.dfb_auto_scale)
            self.chk_dfb_adaptive.setChecked(cfg.dfb_adaptive)
            self.f_dfb_I1.setValue(cfg.dfb_I1_mA)
            self.f_dfb_G1.setValue(cfg.dfb_G1)
            self.f_dfb_I2.setValue(cfg.dfb_I2_mA)
            self.f_dfb_G2.setValue(cfg.dfb_G2)

            self.chk_hom_align.setChecked(cfg.hom_auto_align)
            self.hom_table.setRowCount(0)
            for m in cfg.homs:
                self._append_hom_row(m.Q, m.R, m.freq_Hz)
            if hasattr(self, "chk_map_ff_iter"):
                self.chk_map_ff_iter.blockSignals(True)
                self.chk_map_ff_iter.setChecked(cfg.use_form_factor_iter)
                self.chk_map_ff_iter.blockSignals(False)
        finally:
            self._updating_ui_batch = False
        self._update_estimated_walltime_ui()
        self._update_natural_bunch_length_label()
        self._sync_dmode_prescreen_availability()

    def _update_natural_bunch_length_label(self) -> None:
        try:
            c_val = self.f_C.value()
            e0 = self.f_E0_GeV.value() * 1e9
            u0 = self.f_U0_keV.value() * 1e3
            v_mc = self.f_V_mc_MV.value() * 1e6
            h_val = self.f_h.value()
            ac = self.f_alpha_c.value()
            se0 = self.f_sigma_e0.value()
            qs = math.sqrt(max(1e-12, (ac / (2.0 * math.pi)) * h_val * math.sqrt(max(0.0, (v_mc / e0) ** 2 - (u0 / e0) ** 2))))
            t0 = c_val / 299792458.0
            st0_ps = (ac * se0 / (2.0 * math.pi * qs)) * t0 * 1e12
            self.lbl_sigma_t0.setText(f"{st0_ps:.2f} ps")
        except Exception:
            pass

    def _on_ff_iter_toggled(self, val: bool, source: str) -> None:
        if self._updating_ui_batch:
            return
        self._updating_ui_batch = True
        try:
            if source == "map":
                self.chk_ff_iter.setChecked(val)
            elif hasattr(self, "chk_map_ff_iter"):
                self.chk_map_ff_iter.setChecked(val)
        finally:
            self._updating_ui_batch = False
        self._redraw_stability_map()
        self._refresh_nyquist_plot()

    def _on_optimize_hc_detuning(self) -> None:
        """Exact port of MATLAB `OptimizeHCDetuningButtonPushed` (STABLE_App.m lines 1897-1981)."""
        try:
            cfg = self._collect_machine_config()
            target_xi = 0.98
            fp_opt = calc_bunch_form_factor_eq(
                "detune_for_xi",
                cfg.I0,
                cfg.V_mc,
                cfg.U0,
                cfg.E0,
                cfg.C,
                cfg.h,
                cfg.alpha_c,
                cfg.sigma_e0,
                cfg.n_hc,
                cfg.R_hc,
                cfg.Q_hc,
                target_xi,
                cfg.use_form_factor_iter,
            )
            opt_df_hz = float(fp_opt["fre_shift_hc"]) if float(fp_opt["fre_shift_hc"]) > 0 else 35e3
            opt_df_khz = opt_df_hz * 1e-3
            opt_xi = float(fp_opt["xi_eff"])
            st0 = float(fp_opt.get("sigma_t0", 0.0))
            st_eq = float(fp_opt.get("sigma_t_eq", 0.0))
            opt_stretch = (st_eq / st0) if (cfg.use_form_factor_iter and st0 > 0) else (1.0 + 3.2 * (opt_xi ** 2))

            self.chk_hc_auto.setChecked(False)
            self.f_df_hc_kHz.setValue(opt_df_khz)
            self._set_status(
                f"Optimized: Df_hc = {opt_df_khz:.2f} kHz (xi={opt_xi:.2f}, STABLE, Bunch Stretch: {opt_stretch:.2f}x)"
            )
            QMessageBox.information(
                self,
                "HC Detuning Optimized",
                "Physics Stability-Constrained Optimization Succeeded!\n\n"
                f"Target Current: {cfg.I0 * 1e3:.1f} mA\n"
                f"Optimal Detuning: {opt_df_khz:.2f} kHz ({opt_df_hz:.0f} Hz)\n"
                f"Flat Potential Parameter xi: {opt_xi:.3f}\n"
                f"Predicted Bunch Lengthening: {opt_stretch:.2f}x\n"
                "Stability Status: STABLE",
            )
        except Exception as exc:
            QMessageBox.warning(self, "Optimization Error", str(exc))

    def _on_align_homs_button(self) -> None:
        """Exact port of MATLAB `AlignHOMButtonPushed` (STABLE_App.m lines 2201-2222)."""
        cfg = self._collect_machine_config()
        if not cfg.homs:
            return
        aligned = align_homs_to_worst_case_sideband(cfg)
        self.hom_table.setRowCount(0)
        for m in aligned:
            self._append_hom_row(m.Q, m.R, m.freq_Hz)
        self._set_status("HOM frequencies aligned to worst-case sideband (mbtrack2 standard).")

    def _on_load_config_file(self) -> None:
        """Exact port of MATLAB `LoadConfigButtonPushed` (.mat or .m machine file)."""
        fpath, _ = QFileDialog.getOpenFileName(
            self,
            "Load STABLE Configuration or Machine File",
            "",
            "STABLE Config or Machine (*.mat *.m);;MAT Files (*.mat);;MATLAB Machine Script (*.m)",
        )
        if not fpath:
            return
        try:
            m_cfg, s_cfg = load_config_from_mat_or_m(
                fpath, base_cfg=self._collect_machine_config(), base_scan=self._collect_scan_config()
            )
            self._apply_machine_config_to_ui(m_cfg)
            self.f_scan_df_min.setValue(s_cfg.Df_min_kHz)
            self.f_scan_df_max.setValue(s_cfg.Df_max_kHz)
            self.f_scan_df_step.setValue(s_cfg.Df_step_kHz)
            self.f_scan_i0_min.setValue(s_cfg.I0_min_mA)
            self.f_scan_i0_max.setValue(s_cfg.I0_max_mA)
            self.f_scan_i0_step.setValue(s_cfg.I0_step_mA)
            self._redraw_stability_map()
            self._refresh_nyquist_plot()
            self._set_status(f"Loaded configuration: {Path(fpath).name}")
        except Exception as exc:
            QMessageBox.warning(self, "Load Config Error", str(exc))

    def _on_save_config_file(self) -> None:
        """Exact port of MATLAB `SaveConfigButtonPushed` (.mat config file)."""
        fpath, _ = QFileDialog.getSaveFileName(
            self, "Save STABLE Configuration (.mat)", "stable_config.mat", "MAT Files (*.mat)"
        )
        if not fpath:
            return
        try:
            mat_dict = machine_config_to_matlab_struct(self._collect_machine_config(), self._collect_scan_config())
            sio.savemat(fpath, mat_dict)
            self._set_status(f"Saved configuration: {Path(fpath).name}")
        except Exception as exc:
            QMessageBox.warning(self, "Save Config Error", str(exc))

    def _collect_scan_config(self) -> ScanConfig:
        return ScanConfig(
            Df_min_kHz=self.f_scan_df_min.value(),
            Df_max_kHz=self.f_scan_df_max.value(),
            Df_step_kHz=self.f_scan_df_step.value(),
            I0_min_mA=self.f_scan_i0_min.value(),
            I0_max_mA=self.f_scan_i0_max.value(),
            I0_step_mA=self.f_scan_i0_step.value(),
            enable_dmode_prescreen=self.chk_dmode_prescreen.isChecked(),
            num_gpus=self.f_scan_gpus.value(),
            workers_per_gpu=self.f_workers_per_gpu.value(),
            execution_mode=self.f_exec_mode.currentText(),
        )

    def _update_estimated_walltime_ui(self) -> None:
        if self._updating_ui_batch:
            return
        n_g = self.f_scan_gpus.value()
        w_g = self.f_workers_per_gpu.value()
        vram_gb = round(0.35 * w_g, 1)
        self.lbl_gpu_concurrency.setText(
            f"{n_g} A100 × {w_g} Pts/GPU = {n_g * w_g} points in parallel (~{vram_gb} GB / 80 GB VRAM per GPU)"
        )
        auto_on = self.chk_auto_walltime.isChecked()
        self.f_walltime_hours.setEnabled(not auto_on)
        if auto_on:
            try:
                m_cfg = self._collect_machine_config()
                s_cfg = self._collect_scan_config()
                tmp_prof = TGCCConnectionProfile(auto_walltime=True)
                est_sec = tmp_prof.compute_effective_walltime_seconds(m_cfg, s_cfg, is_scan=True)
                self.f_walltime_hours.blockSignals(True)
                self.f_walltime_hours.setValue(round(est_sec / 3600.0, 1))
                self.f_walltime_hours.blockSignals(False)
            except Exception:
                pass

    def _collect_tgcc_profile(self) -> TGCCConnectionProfile:
        return TGCCConnectionProfile(
            profile_name="TGCC Irene A100",
            mode=self.f_conn_mode.currentText(),
            ssh_host=self.f_ssh_host.text().strip(),
            ssh_port=self.f_ssh_port.value(),
            ssh_username=self.f_ssh_user.text().strip(),
            ccc_project=self.f_ccc_proj.text().strip(),
            ccc_partition=self.f_ccc_part.text().strip(),
            ccc_qos=self.f_ccc_qos.currentText(),
            remote_work_dir=self.f_remote_dir.text().strip(),
            python_env_cmd=self.f_python_cmd.text().strip(),
            walltime_seconds=int(round(self.f_walltime_hours.value() * 3600.0)),
            auto_walltime=self.chk_auto_walltime.isChecked(),
        )

    def _load_user_profile_and_presets(self) -> None:
        prof, secret = self.vault.load_tgcc_profile(self.session)
        if prof.ssh_host in ("", "irene-fr.ccc.cea.fr"):
            prof.ssh_host = "topaze.ccc.cea.fr"
        if prof.remote_work_dir in ("", "$CCCSCRATCHDIR/STABLE_Py_TGCC"):
            prof.remote_work_dir = "/ccc/work/cont003/soleil/zhaolu/STABLE_Py_TGCC"
        if prof.ccc_project in ("", "gen1234", "soleil"):
            prof.ccc_project = "soldai"
        if not prof.ssh_username:
            prof.ssh_username = "zhaolu"

        self.f_conn_mode.setCurrentText(prof.mode)
        self.f_ssh_host.setText(prof.ssh_host)
        self.f_ssh_port.setValue(prof.ssh_port)
        self.f_ssh_user.setText(prof.ssh_username)
        self.f_ssh_pwd.setText(secret)
        self.f_ccc_proj.setText(prof.ccc_project)
        self.f_ccc_part.setText(prof.ccc_partition)
        self.f_ccc_qos.setCurrentText(getattr(prof, "ccc_qos", "normal") or "normal")
        self.f_remote_dir.setText(prof.remote_work_dir)
        self.f_python_cmd.setText(prof.python_env_cmd)
        self.chk_auto_walltime.setChecked(prof.auto_walltime)
        self.f_walltime_hours.setValue(max(0.5, min(24.0, round(prof.walltime_seconds / 3600.0, 1))))
        self._update_estimated_walltime_ui()

        presets = self.vault.list_user_presets(self.session)
        self.preset_combo.blockSignals(True)
        self.preset_combo.clear()
        self.preset_combo.addItems(presets)
        self.preset_combo.blockSignals(False)

        last_state = None
        try:
            last_state = self.vault.load_last_session_state(self.session)
        except Exception:
            last_state = None

        if last_state is not None:
            m_cfg, s_cfg = last_state
            self._apply_scan_config_to_ui(s_cfg)
            self._apply_machine_config_to_ui(m_cfg)
            if m_cfg.preset_name and self.preset_combo.findText(m_cfg.preset_name) >= 0:
                self.preset_combo.blockSignals(True)
                self.preset_combo.setCurrentText(m_cfg.preset_name)
                self.preset_combo.blockSignals(False)
            self._append_log("[Session] Restored settings from previous session.")
        elif presets:
            self._on_preset_selected(presets[0])

    def _apply_scan_config_to_ui(self, s_cfg: ScanConfig) -> None:
        self.f_scan_df_min.setValue(s_cfg.Df_min_kHz)
        self.f_scan_df_max.setValue(s_cfg.Df_max_kHz)
        self.f_scan_df_step.setValue(s_cfg.Df_step_kHz)
        self.f_scan_i0_min.setValue(s_cfg.I0_min_mA)
        self.f_scan_i0_max.setValue(s_cfg.I0_max_mA)
        self.f_scan_i0_step.setValue(s_cfg.I0_step_mA)
        self.f_workers_per_gpu.setValue(s_cfg.workers_per_gpu if s_cfg.workers_per_gpu > 5 else 16)
        self.chk_dmode_prescreen.setChecked(bool(s_cfg.enable_dmode_prescreen))
        self._update_estimated_walltime_ui()

    def _sync_dmode_prescreen_availability(self) -> None:
        """Disable D-mode Robinson pre-screening for normal-conducting cavities (Q_hc < 1e5, e.g. SOLEIL II)."""
        q_hc = float(self.f_Q_hc.value()) if hasattr(self, "f_Q_hc") else 31000.0
        is_sc = q_hc >= 1.0e5
        if not is_sc:
            self.chk_dmode_prescreen.setChecked(False)
            self.chk_dmode_prescreen.setEnabled(False)
            self.chk_dmode_prescreen.setText("D-Mode Robinson Pre-Screening (SC Q_hc >= 1e5 only; Tracking for NC)")
        else:
            self.chk_dmode_prescreen.setEnabled(True)
            self.chk_dmode_prescreen.setText("Enable D-Mode Robinson Pre-Screening (T. He PRAB 2024, SC)")

    def _on_apply_feedback_preset(self) -> None:
        """Apply quick feedback preset for xi=1 stability island tuning and refresh Nyquist/Bode."""
        idx = self.fb_preset_combo.currentIndex()
        self._updating_ui_batch = True
        try:
            if idx == 0:
                # SOLEIL II Default
                self.chk_pi_enable.setChecked(True)
                self.f_pi_kp.setValue(0.8)
                self.f_pi_ki.setValue(1.0e-4)
                self.f_pi_m.setValue(80)
                self.f_pi_d.setValue(800)
                self.chk_dfb_enable.setChecked(True)
                self.f_dfb_gain.setValue(0.20)
                self.f_dfb_phase.setValue(0.0)
                self.f_dfb_delay.setValue(704)
                self.f_dfb_m.setValue(208)
                self.f_dfb_every.setValue(208)
                self.chk_dfb_autoscale.setChecked(True)
                self.chk_dfb_adaptive.setChecked(False)
                self.f_dfb_I1.setValue(300.0)
                self.f_dfb_G1.setValue(0.25)
                self.f_dfb_I2.setValue(500.0)
                self.f_dfb_G2.setValue(0.20)
            elif idx == 1:
                # mbtrack2 Weak DFB (G=0.15)
                self.chk_pi_enable.setChecked(True)
                self.f_pi_kp.setValue(0.8)
                self.f_pi_ki.setValue(1.0e-4)
                self.chk_dfb_enable.setChecked(True)
                self.chk_dfb_autoscale.setChecked(False)
                self.chk_dfb_adaptive.setChecked(False)
                self.f_dfb_gain.setValue(0.15)
                self.f_dfb_phase.setValue(0.0)
            elif idx == 2:
                # mbtrack2 Moderate DFB (G=0.25)
                self.chk_pi_enable.setChecked(True)
                self.f_pi_kp.setValue(1.0)
                self.f_pi_ki.setValue(1.0e-4)
                self.chk_dfb_enable.setChecked(True)
                self.chk_dfb_autoscale.setChecked(False)
                self.chk_dfb_adaptive.setChecked(False)
                self.f_dfb_gain.setValue(0.25)
                self.f_dfb_phase.setValue(0.0)
            elif idx == 3:
                # PI Only
                self.chk_pi_enable.setChecked(True)
                self.f_pi_kp.setValue(0.8)
                self.f_pi_ki.setValue(1.0e-4)
                self.chk_dfb_enable.setChecked(False)
            elif idx == 4:
                # Open Loop
                self.chk_pi_enable.setChecked(False)
                self.chk_dfb_enable.setChecked(False)
        finally:
            self._updating_ui_batch = False
        self._refresh_nyquist_plot()
        self._set_status(f"Applied feedback preset: {self.fb_preset_combo.currentText()}")
        self._append_log(f"[Feedback] Applied '{self.fb_preset_combo.currentText()}'")

    def _on_save_to_soleil_default(self) -> None:
        """Save current full MachineConfig + ScanConfig into 'SOLEIL II (Default)' preset."""
        target = "SOLEIL II (Default)"
        m_cfg = self._collect_machine_config()
        s_cfg = self._collect_scan_config()
        self.vault.save_user_preset(self.session, target, m_cfg, s_cfg)
        self.vault.save_last_session_state(self.session, m_cfg, s_cfg)
        presets = self.vault.list_user_presets(self.session)
        self.preset_combo.blockSignals(True)
        self.preset_combo.clear()
        self.preset_combo.addItems(presets)
        self.preset_combo.setCurrentText(target)
        self.preset_combo.blockSignals(False)
        self._set_status(f"Saved current settings to '{target}'.")
        self._append_log(f"[Preset] Saved current settings into '{target}'.")

    def _prompt_save_settings_after_tracking(self, source_label: str) -> None:
        """Ask the user after tracking completes whether to save current settings to SOLEIL II Default or a new preset."""
        box = QMessageBox(self)
        box.setWindowTitle(f"{source_label} Finished — Save Settings?")
        box.setText(
            f"{source_label} completed.\n"
            "Would you like to save the current machine & feedback settings?"
        )
        btn_soleil = box.addButton("Save to SOLEIL II (Default)", QMessageBox.ButtonRole.AcceptRole)
        btn_new = box.addButton("Save as New Preset...", QMessageBox.ButtonRole.ActionRole)
        box.addButton("Keep for Session Only", QMessageBox.ButtonRole.RejectRole)
        box.exec()
        clicked = box.clickedButton()
        if clicked == btn_soleil:
            self._on_save_to_soleil_default()
        elif clicked == btn_new:
            self._on_save_preset()

    def _on_preset_selected(self, name: str) -> None:
        if not name:
            return
        try:
            m_cfg, s_cfg = self.vault.load_user_preset(self.session, name)
            self._apply_machine_config_to_ui(m_cfg)
            self.f_scan_df_min.setValue(s_cfg.Df_min_kHz)
            self.f_scan_df_max.setValue(s_cfg.Df_max_kHz)
            self.f_scan_df_step.setValue(s_cfg.Df_step_kHz)
            self.f_scan_i0_min.setValue(s_cfg.I0_min_mA)
            self.f_scan_i0_max.setValue(s_cfg.I0_max_mA)
            self.f_scan_i0_step.setValue(s_cfg.I0_step_mA)
            if s_cfg.workers_per_gpu > 5:
                self.f_workers_per_gpu.setValue(s_cfg.workers_per_gpu)
            else:
                self.f_workers_per_gpu.setValue(16)
            self._update_estimated_walltime_ui()
            self._append_log(f"[Preset] Loaded user preset '{name}'")
        except Exception as exc:
            self._append_log(f"[Preset Error] {exc}")

    def _on_save_preset(self) -> None:
        name, ok = QInputDialog.getText(
            self, "Save User Preset", "Preset Name:", text=self.preset_combo.currentText()
        )
        if ok and name.strip():
            m_cfg = self._collect_machine_config()
            s_cfg = self._collect_scan_config()
            self.vault.save_user_preset(self.session, name.strip(), m_cfg, s_cfg)
            presets = self.vault.list_user_presets(self.session)
            self.preset_combo.blockSignals(True)
            self.preset_combo.clear()
            self.preset_combo.addItems(presets)
            self.preset_combo.setCurrentText(name.strip())
            self.preset_combo.blockSignals(False)
            self._append_log(f"[Preset] Saved preset '{name.strip()}' for user {self.session.username}")

    def _on_delete_preset(self) -> None:
        cur = self.preset_combo.currentText()
        if not cur:
            return
        self.vault.delete_user_preset(self.session, cur)
        self._load_user_profile_and_presets()

    def _on_save_tgcc_profile(self) -> None:
        prof = self._collect_tgcc_profile()
        secret = self.f_ssh_pwd.text()
        self.vault.save_tgcc_profile(self.session, prof, ssh_secret=secret)
        QMessageBox.information(
            self,
            "Encrypted Profile Saved",
            f"TGCC SSH profile & password encrypted and saved for user '{self.session.username}'.",
        )

    # =========================================================================
    # Backend Connection & TGCC Deployment (Asynchronous + Animated Spinner)
    # =========================================================================
    def _set_top_busy(self, busy: bool, badge_text: Optional[str] = None) -> None:
        if busy:
            self.top_spinner.setVisible(True)
            self.top_spinner.start()
            if badge_text:
                self.conn_badge.setText(badge_text)
                self.conn_badge.setStyleSheet(
                    "background: #e0f2fe; color: #0369a1; padding: 4px 8px; border-radius: 4px; font-weight: bold;"
                )
        else:
            self.top_spinner.stop()
            self.top_spinner.setVisible(False)

    def _on_async_task_step(self, msg: str, pct: int) -> None:
        self._set_status(msg)
        if self._busy_dialog is not None:
            self._busy_dialog.update_status(msg, pct)

    def _connect_backend(self, silent: bool = False) -> None:
        prof = self._collect_tgcc_profile()
        secret = self.f_ssh_pwd.text()
        otp = self.f_ssh_otp.text().strip()
        if silent and prof.mode == "ssh_tgcc" and not secret:
            prof = TGCCConnectionProfile(mode="local")

        target_label = f"{prof.ssh_username}@{prof.ssh_host}" if prof.mode == "ssh_tgcc" else "Local Physics Daemon"
        self._set_top_busy(True, f"Connecting to {target_label}...")

        if not silent:
            self._busy_dialog = AsyncTaskProgressDialog(
                title=f"Connecting to {target_label}",
                initial_text=f"Opening connection to {target_label}...",
                parent=self,
            )
            self._busy_dialog.show()

        def _worker():
            try:
                st = self.connector.connect(
                    prof,
                    ssh_secret=secret,
                    otp_code=otp,
                    log_fn=lambda m: self.signals.log_message.emit(m),
                    step_fn=lambda s, p: self.signals.async_task_step.emit(s, p),
                )
                info = ""
                if not silent and prof.mode == "ssh_tgcc":
                    self.signals.async_task_step.emit("[4/4] Querying Topaze A100 partition & ccc_mpp queue...", 92)
                    info = self.connector.query_tgcc_cluster_info()
                self.signals.connect_finished.emit(st, info, silent, prof)
            except Exception as exc:
                self.signals.connect_failed.emit(str(exc), silent)

        threading.Thread(target=_worker, daemon=True).start()

    def _on_connect_finished(self, st: Dict[str, Any], info: str, silent: bool, prof: TGCCConnectionProfile) -> None:
        if self._busy_dialog is not None:
            self._busy_dialog.finish_and_close()
            self._busy_dialog = None
        self._set_top_busy(False)

        gpus = st.get("gpus", [])
        hname = st.get("hostname", "server")
        if prof.mode == "ssh_tgcc":
            badge_str = f"Gateway: {hname} -> Queue: {prof.ccc_partition} (NVIDIA A100-80GB, Account: {prof.ccc_project})"
        else:
            gpu_desc = f"{len(gpus)}x {gpus[0]['name']}" if gpus else "CPU"
            badge_str = f"Local Server: {hname} [{gpu_desc}]"
        self.conn_badge.setText(badge_str)
        self.conn_badge.setStyleSheet("background: #dcfce7; color: #166534; padding: 4px 8px; border-radius: 4px; font-weight: bold;")
        self._set_status(f"Connected: {badge_str}")

        if info:
            self._append_log(info)
        if not silent:
            QMessageBox.information(
                self,
                "Connected to Server",
                f"Successfully connected to {hname}!\n"
                f"Work Dir: {self.connector._remote_work_dir}\n"
                f"Dispatch Mode: {badge_str}",
            )

    def _on_connect_failed(self, err_msg: str, silent: bool) -> None:
        if self._busy_dialog is not None:
            self._busy_dialog.finish_and_close()
            self._busy_dialog = None
        self._set_top_busy(False)

        self.conn_badge.setText("Server: Disconnected")
        self.conn_badge.setStyleSheet("background: #fee2e2; color: #991b1b; padding: 4px 8px; border-radius: 4px; font-weight: bold;")
        self._set_status(f"Connection failed: {err_msg}")
        if not silent:
            QMessageBox.warning(self, "Connection Error", err_msg)

    def _on_deploy_to_tgcc(self) -> None:
        prof = self._collect_tgcc_profile()
        self._set_top_busy(True, f"Deploying code to {prof.ssh_host}...")
        self._busy_dialog = AsyncTaskProgressDialog(
            title="Deploying STABLE_Py_TGCC to Topaze",
            initial_text=f"Packing and streaming Python package to {prof.ssh_host}:{prof.remote_work_dir}...",
            parent=self,
        )
        self._busy_dialog.show()

        def _worker():
            try:
                msg = self.connector.deploy_package_to_tgcc(
                    prof,
                    log_fn=lambda m: self.signals.log_message.emit(m),
                    force=True,
                )
                self.signals.deploy_finished.emit(msg)
            except Exception as exc:
                self.signals.deploy_failed.emit(str(exc))

        threading.Thread(target=_worker, daemon=True).start()

    def _on_deploy_finished(self, msg: str) -> None:
        if self._busy_dialog is not None:
            self._busy_dialog.finish_and_close()
            self._busy_dialog = None
        self._set_top_busy(False)
        self._set_status(msg)
        QMessageBox.information(self, "TGCC Code Deployment Complete", msg)

    def _on_deploy_failed(self, err: str) -> None:
        if self._busy_dialog is not None:
            self._busy_dialog.finish_and_close()
            self._busy_dialog = None
        self._set_top_busy(False)
        self._set_status(f"Deploy failed: {err}")
        QMessageBox.warning(self, "Deploy Failed", err)

    def _on_clean_pending_queue(self) -> None:
        if self.connector.active_mode != "ssh_tgcc":
            QMessageBox.information(self, "Queue Cleanup", "Connect to TGCC via SSH first.")
            return
        self._set_status("Cleaning leftover STABLE_* jobs on Topaze...")

        def _worker():
            try:
                out = self.connector.clean_stale_stable_jobs(log_fn=lambda m: self.signals.log_message.emit(m))
                info = self.connector.query_tgcc_cluster_info()
                self.signals.log_message.emit(f"[Queue Cleanup Result]\n{out}\n{info}")
                self.signals.status_text.emit("Cleaned pending STABLE_* jobs on Topaze (see Console Tab 4).")
            except Exception as exc:
                self.signals.job_error.emit(str(exc))

        threading.Thread(target=_worker, daemon=True).start()

    # =========================================================================
    # Execution: Single Point Tracking & 2D Multi-A100 Scan
    # =========================================================================
    def _start_single_tracking(self) -> None:
        cfg = self._collect_machine_config()
        prof = self._collect_tgcc_profile()
        exec_mode = self.f_exec_mode.currentText()
        job_id = f"track_{int(time.time())}"
        self._current_job_id = job_id
        self.progress_bar.setValue(0)
        self.right_tabs.setCurrentIndex(0)
        self._set_top_busy(True)
        self._set_status(f"Running single tracking {job_id} ({exec_mode})...")

        def worker():
            self.connector.stream_single_tracking(
                cfg=cfg,
                job_id=job_id,
                on_frame=lambda fr: self.signals.telemetry_received.emit(fr),
                on_error=lambda err: self.signals.job_error.emit(err),
                execution_mode=exec_mode,
                profile=prof,
                log_fn=lambda m: self.signals.log_message.emit(m),
            )

        threading.Thread(target=worker, daemon=True).start()

    def _start_2d_scan(self) -> None:
        cfg = self._collect_machine_config()
        scan_cfg = self._collect_scan_config()
        prof = self._collect_tgcc_profile()
        job_id = f"scan_{int(time.time())}"
        self._current_job_id = job_id

        # Initialize 2D scan matrices matching MATLAB STABLE_App.m
        self._df_grid_kHz = np.array(scan_cfg.get_df_grid_kHz(), dtype=np.float64)
        self._i0_grid_mA = np.array(scan_cfg.get_i0_grid_mA(), dtype=np.float64)
        n_i0 = len(self._i0_grid_mA)
        n_df = len(self._df_grid_kHz)
        self._stability_code_mat = np.full((n_i0, n_df), np.nan, dtype=np.float64)
        self._growth_rate_mat = np.full((n_i0, n_df), np.nan, dtype=np.float64)
        self._bunch_len_mat = np.full((n_i0, n_df), np.nan, dtype=np.float64)
        self._jitter_mat = np.full((n_i0, n_df), np.nan, dtype=np.float64)
        self._show_critical_boundary_dots = False
        self._active_scan_marker = None
        self._active_scan_markers.clear()
        self._completed_count = 0
        self._total_scan_points = n_i0 * n_df

        self.progress_bar.setValue(0)
        self.right_tabs.setCurrentIndex(1)
        self._redraw_stability_map()
        self._set_top_busy(True)
        self._set_status(
            f"Submitting & Running 2D Scan ({n_i0}x{n_df} = {self._total_scan_points} points, {scan_cfg.num_gpus} GPUs)..."
        )

        def worker():
            self.connector.stream_2d_scan(
                cfg=cfg,
                scan_cfg=scan_cfg,
                profile=prof,
                job_id=job_id,
                on_point=lambda pt: self.signals.scan_point_received.emit(pt),
                on_progress=lambda pr: self.signals.scan_progress_received.emit(pr),
                on_done=lambda n: self.signals.scan_completed.emit(n),
                on_error=lambda err: self.signals.job_error.emit(err),
                log_fn=lambda m: self.signals.log_message.emit(m),
            )

        threading.Thread(target=worker, daemon=True).start()

    def _stop_active_job(self) -> None:
        self._set_top_busy(False)
        self._set_status("Stopping active job(s)...")
        jid = self._current_job_id or "*"
        threading.Thread(target=lambda: self.connector.stop_job(jid), daemon=True).start()

    def _auto_load_latest_scan_on_startup(self) -> None:
        """Port of MATLAB `updateStabilityMapView` auto-discovery on startup."""
        try:
            latest = find_latest_default_scan_summary()
            if latest is not None and latest.exists():
                self._load_scan_summary_path(latest, restore_config=False, switch_tab=False)
            else:
                self._redraw_stability_map()
        except Exception:
            self._redraw_stability_map()

    def _load_scan_summary_path(
        self,
        fpath: str | Path,
        restore_config: bool = True,
        switch_tab: bool = True,
    ) -> None:
        """
        Load `scan_summary.mat` and (when `restore_config=True`) restore its exact `MachineConfig`
        and `ScanConfig` to the UI (`Q_hc`, `R_hc`, `V_mc`, `PI`, `DFB`, `HOMs`) just like MATLAB
        `loadScanDataForMap` + `restoreScanConfigToUI` (STABLE_App.m lines 1670-1698 & 2588-2720).
        """
        p = Path(fpath)
        s = sio.loadmat(str(p), squeeze_me=True, struct_as_record=False)
        if "Df_scan_kHz" not in s or "I0_scan_mA" not in s or "stability_code" not in s:
            raise ValueError("scan_summary.mat is missing required fields (Df_scan_kHz, I0_scan_mA, stability_code).")

        if restore_config:
            m_cfg, s_cfg = restore_scan_config_from_mat(
                s, p, base_cfg=MachineConfig.soleil2_default(), base_scan=self._collect_scan_config()
            )
            self._apply_machine_config_to_ui(m_cfg)
            self.f_scan_df_min.setValue(s_cfg.Df_min_kHz)
            self.f_scan_df_max.setValue(s_cfg.Df_max_kHz)
            self.f_scan_df_step.setValue(s_cfg.Df_step_kHz)
            self.f_scan_i0_min.setValue(s_cfg.I0_min_mA)
            self.f_scan_i0_max.setValue(s_cfg.I0_max_mA)
            self.f_scan_i0_step.setValue(s_cfg.I0_step_mA)
            self._refresh_nyquist_plot()

        self._df_grid_kHz = np.atleast_1d(np.asarray(s["Df_scan_kHz"], dtype=np.float64)).ravel()
        self._i0_grid_mA = np.atleast_1d(np.asarray(s["I0_scan_mA"], dtype=np.float64)).ravel()
        sc_raw = np.asarray(s["stability_code"], dtype=np.float64)
        if sc_raw.ndim == 1:
            sc_raw = sc_raw.reshape(len(self._i0_grid_mA), len(self._df_grid_kHz))
        self._stability_code_mat = sc_raw

        if "growth_rate" in s and np.asarray(s["growth_rate"]).size == self._stability_code_mat.size:
            gr_raw = np.asarray(s["growth_rate"], dtype=np.float64).reshape(self._stability_code_mat.shape)
        else:
            gr_raw = self._stability_code_mat * 50.0
        self._growth_rate_mat = gr_raw
        self._bunch_len_mat = np.full_like(self._stability_code_mat, np.nan)
        self._jitter_mat = np.full_like(self._stability_code_mat, np.nan)
        self._show_critical_boundary_dots = False
        self._active_scan_marker = None
        self._active_scan_markers.clear()
        if switch_tab:
            self.right_tabs.setCurrentIndex(1)
        self._redraw_stability_map()
        self._set_status(f"Map & Settings loaded: {p.parent.name}/{p.name}")
        self._append_log(f"[Loaded Scan Summary] {p}")

    def _on_load_matlab_scan_summary(self) -> None:
        default_dir = str(Path(__file__).resolve().parents[3] / "scan_data")
        fpath, _ = QFileDialog.getOpenFileName(
            self, "Load MATLAB scan_summary.mat", default_dir, "MAT Files (*.mat)"
        )
        if not fpath:
            return
        try:
            self._load_scan_summary_path(fpath, restore_config=True, switch_tab=True)
        except Exception as exc:
            QMessageBox.warning(self, "Load .mat Error", str(exc))

    def _on_load_scan_folder(self) -> None:
        """Exact port of MATLAB `loadScanDataForMap` folder picker (`uigetdir`)."""
        default_dir = str(Path(__file__).resolve().parents[3] / "scan_data")
        folder = QFileDialog.getExistingDirectory(self, "Select Scan Data Folder (Scan_*)", default_dir)
        if not folder:
            return
        sum_file = Path(folder) / "scan_summary.mat"
        if not sum_file.exists():
            QMessageBox.warning(self, "File Not Found", f"No scan_summary.mat found in:\n{folder}")
            return
        try:
            self._load_scan_summary_path(sum_file, restore_config=True, switch_tab=True)
        except Exception as exc:
            QMessageBox.warning(self, "Load Folder Error", str(exc))

    def _on_save_scan_summary_mat(self) -> None:
        """Save current 2D scan map and UI settings to a MATLAB-compatible `scan_summary.mat`."""
        if self._df_grid_kHz.size == 0 or self._i0_grid_mA.size == 0:
            QMessageBox.information(self, "No Scan Data", "Run a 2D scan or load scan data first.")
            return
        fpath, _ = QFileDialog.getSaveFileName(
            self, "Save MATLAB scan_summary.mat", "scan_summary.mat", "MAT Files (*.mat)"
        )
        if not fpath:
            return
        try:
            cfg = self._collect_machine_config()
            s_cfg = self._collect_scan_config()
            cfg_struct = machine_config_to_matlab_struct(cfg, s_cfg)
            dfb_gains = np.array([cfg.get_dfb_gain_for_current(i0 * 1e-3) for i0 in self._i0_grid_mA], dtype=np.float64)
            sio.savemat(
                fpath,
                {
                    "Df_scan_kHz": self._df_grid_kHz,
                    "I0_scan_mA": self._i0_grid_mA,
                    "DFB_gains_mA": dfb_gains,
                    "stability_code": self._stability_code_mat,
                    "growth_rate": self._growth_rate_mat,
                    "config": cfg_struct,
                    "scanFolderName": Path(fpath).parent.name,
                },
            )
            self._set_status(f"Saved MATLAB scan summary: {fpath}")
        except Exception as exc:
            QMessageBox.warning(self, "Save Error", str(exc))

    def _on_fit_safety_boundary(self) -> None:
        """Exact port of MATLAB `fitAndPlotSafetyBoundary` (STABLE_App.m lines 1984-2069)."""
        if self._stability_code_mat.size == 0 or np.all(np.isnan(self._stability_code_mat)):
            QMessageBox.information(self, "No Scan Data", "No scan data loaded to extract safety boundary.")
            return
        fit_res = fit_safety_boundary(self._df_grid_kHz, self._i0_grid_mA, self._stability_code_mat)
        if fit_res is None:
            self._set_status("Notice: No clear boundary transition in current scan data.")
            QMessageBox.information(self, "Safety Boundary", "No clear unstable -> stable transition found across rows.")
            return
        self._show_critical_boundary_dots = True
        self._redraw_stability_map()
        self._set_status(fit_res["status_msg"])
        QMessageBox.information(self, "Safety Boundary Extracted", fit_res["dialog_msg"])

    # =========================================================================
    # Live Plotting & Telemetry Handlers
    # =========================================================================
    def _on_live_view_mode_changed(self) -> None:
        self._line_hc_amp = None
        if self._last_telemetry_frame is not None:
            fr = self._last_telemetry_frame
            self._update_four_panel_waveforms(
                fr.turn, fr.hc_amp_kV, fr.hc_phase_deg, fr.mc_amp_kV, fr.mc_phase_deg, force=True
            )

    def _on_telemetry_frame(self, frame: TelemetryFrame) -> None:
        self._last_telemetry_frame = frame
        pct = int(round(100.0 * frame.turn / max(1, frame.total_turns)))
        self.progress_bar.setValue(pct)
        self._set_status(
            f"[{frame.device_name}] Turn {frame.turn}/{frame.total_turns} ({pct}%) | "
            f"xi_eff={frame.xi_eff:.3f}, psi_HC={frame.psi_hc_deg:.1f}° | "
            f"Elapsed: {frame.elapsed_sec:.1f}s, ETA: {frame.eta_sec:.1f}s | {frame.status_label}"
        )
        self._update_four_panel_waveforms(
            frame.turn, frame.hc_amp_kV, frame.hc_phase_deg, frame.mc_amp_kV, frame.mc_phase_deg
        )

        if frame.event_type == "tracking_done":
            self._set_top_busy(False)
            self._update_four_panel_waveforms(
                frame.turn, frame.hc_amp_kV, frame.hc_phase_deg, frame.mc_amp_kV, frame.mc_phase_deg, force=True
            )
            self._append_log(
                f"[Tracking Done] {frame.status_label} (code={frame.stability_code}) | "
                f"gr={frame.growth_rate_s1:+.2f} s^-1 | f_dom={frame.dominant_freq_Hz:.1f} Hz | "
                f"mu*={frame.dominant_mode_mu}"
            )
            self.ax_evo.clear()
            self.ax_modes.clear()
            if frame.turns_axis:
                self.ax_evo.plot(frame.turns_axis, frame.bunch_length_ps_history, label="Bunch Length [ps]", lw=1.5)
                self.ax_evo.plot(frame.turns_axis, frame.centroid_jitter_ps_history, label="Centroid Jitter [ps]", lw=1.5)
                self.ax_evo.set_xlabel("Turn")
                self.ax_evo.set_ylabel("[ps]")
                self.ax_evo.set_title("Longitudinal Evolution")
                self.ax_evo.legend(fontsize=8)
                self.ax_evo.grid(True, alpha=0.35)

            if frame.mode_spectrum_amp_ps:
                modes = np.arange(len(frame.mode_spectrum_amp_ps))
                self.ax_modes.stem(modes[:60], frame.mode_spectrum_amp_ps[:60], basefmt=" ")
                self.ax_modes.set_xlabel("Coupled-Bunch Mode Index mu")
                self.ax_modes.set_ylabel("Modal Amplitude [ps]")
                self.ax_modes.set_title(f"Spatial Mode Spectrum (Dominant mu* = {frame.dominant_mode_mu})")
                self.ax_modes.grid(True, alpha=0.35)
            self.canvas_diag.draw_idle()
            self._prompt_save_settings_after_tracking("Single Tracking")

    def _update_four_panel_waveforms(
        self,
        turn: int,
        hc_amp: List[float],
        hc_pha: List[float],
        mc_amp: List[float],
        mc_pha: List[float],
        force: bool = False,
    ) -> None:
        now = time.perf_counter()
        if not force and (now - self._last_waveform_draw_ts) < 0.12:
            return
        self._last_waveform_draw_ts = now

        if hasattr(self, "live_view_combo") and self.live_view_combo.currentIndex() == 1:
            self._line_hc_amp = None
            fr = self._last_telemetry_frame
            for ax in (self.ax_hc_amp, self.ax_mc_amp, self.ax_hc_pha, self.ax_mc_pha):
                ax.clear()
            if fr is not None and fr.turns_axis:
                self.ax_hc_amp.plot(fr.turns_axis, fr.centroid_jitter_ps_history, color="#0284c7", lw=1.4)
                self.ax_hc_amp.set_title(r"Centroid Jitter $\langle \tau \rangle_{\mathrm{rms}}$ [ps]")
                self.ax_hc_amp.set_xlabel("Turns")
                self.ax_hc_amp.grid(True, alpha=0.35)

                self.ax_mc_amp.plot(fr.turns_axis, fr.bunch_length_ps_history, color="#dc2626", lw=1.4)
                self.ax_mc_amp.set_title(r"Bunch Length $\sigma_\tau$ [ps]")
                self.ax_mc_amp.set_xlabel("Turns")
                self.ax_mc_amp.grid(True, alpha=0.35)

                self.ax_hc_pha.plot(fr.turns_axis, fr.energy_spread_history, color="#059669", lw=1.4)
                self.ax_hc_pha.set_title(r"Energy Spread $\sigma_\delta$")
                self.ax_hc_pha.set_xlabel("Turns")
                self.ax_hc_pha.grid(True, alpha=0.35)

                if fr.mode_spectrum_amp_ps:
                    modes = np.arange(min(60, len(fr.mode_spectrum_amp_ps)))
                    self.ax_mc_pha.stem(modes, fr.mode_spectrum_amp_ps[: len(modes)], basefmt=" ")
                    self.ax_mc_pha.set_title(f"Coupled-Bunch Modes (mu* = {fr.dominant_mode_mu})")
                    self.ax_mc_pha.set_xlabel("Mode Index mu")
                    self.ax_mc_pha.grid(True, alpha=0.35)
            self.canvas_live.draw_idle()
            return

        if not hc_amp and not mc_amp:
            return

        if self._line_hc_amp is None or len(hc_amp) != len(self._line_hc_amp.get_xdata()):
            for ax in (self.ax_hc_amp, self.ax_mc_amp, self.ax_hc_pha, self.ax_mc_pha):
                ax.clear()
            x_hc = np.arange(len(hc_amp)) if hc_amp else np.arange(1)
            x_mc = np.arange(len(mc_amp)) if mc_amp else np.arange(1)
            (self._line_hc_amp,) = self.ax_hc_amp.plot(x_hc, hc_amp or [0.0], color="#0284c7", lw=1.2)
            self.ax_hc_amp.set_ylabel("|V_HC| [kV]")
            self.ax_hc_amp.grid(True, alpha=0.35)
            (self._line_hc_pha,) = self.ax_hc_pha.plot(x_hc, hc_pha or [0.0], color="#0284c7", lw=1.2)
            self.ax_hc_pha.set_ylabel("Phase [deg]")
            self.ax_hc_pha.set_xlabel("Bucket Index")
            self.ax_hc_pha.grid(True, alpha=0.35)
            (self._line_mc_amp,) = self.ax_mc_amp.plot(x_mc, mc_amp or [0.0], color="#dc2626", lw=1.2)
            self.ax_mc_amp.set_ylabel("|V_MC| [kV]")
            self.ax_mc_amp.grid(True, alpha=0.35)
            (self._line_mc_pha,) = self.ax_mc_pha.plot(x_mc, mc_pha or [0.0], color="#dc2626", lw=1.2)
            self.ax_mc_pha.set_ylabel("Phase [deg]")
            self.ax_mc_pha.set_xlabel("Bucket Index")
            self.ax_mc_pha.grid(True, alpha=0.35)
        else:
            if hc_amp:
                self._line_hc_amp.set_ydata(hc_amp)
                self._line_hc_pha.set_ydata(hc_pha)
                self.ax_hc_amp.relim()
                self.ax_hc_amp.autoscale_view(scalex=False, scaley=True)
                self.ax_hc_pha.relim()
                self.ax_hc_pha.autoscale_view(scalex=False, scaley=True)
            if mc_amp:
                self._line_mc_amp.set_ydata(mc_amp)
                self._line_mc_pha.set_ydata(mc_pha)
                self.ax_mc_amp.relim()
                self.ax_mc_amp.autoscale_view(scalex=False, scaley=True)
                self.ax_mc_pha.relim()
                self.ax_mc_pha.autoscale_view(scalex=False, scaley=True)

        self.ax_hc_amp.set_title(f"HHC Voltage |V| (Turn {turn})")
        self.ax_mc_amp.set_title(f"MC Voltage |V| (Turn {turn})")
        self.canvas_live.draw_idle()

    def _update_active_scan_marker_fast(self, new_marker: Optional[Tuple[float, float]] = None) -> None:
        """Update ONLY the red active-point ring artist(s) for all concurrent points without rebuilding the 2D heatmap."""
        if new_marker is not None:
            self._active_scan_marker = new_marker
        if self._active_marker_artist is None:
            if not self._map_redraw_timer.isActive():
                self._map_redraw_timer.start()
            return
        try:
            use_psi = self.map_xaxis_combo.currentIndex() == 1
            pts_list = list(self._active_scan_markers.values())
            if not pts_list and self._active_scan_marker is not None:
                pts_list = [self._active_scan_marker]
            if not pts_list:
                self._active_marker_artist.set_data([], [])
            else:
                f0_val = 299792458.0 / max(1.0, self.f_C.value())
                f_hc = self.f_n_hc.value() * self.f_h.value() * f0_val
                q_hc = self.f_Q_hc.value()
                xs = [
                    (math.degrees(math.atan(2.0 * q_hc * (df_v * 1e3) / f_hc)) if use_psi else df_v)
                    for df_v, _ in pts_list
                ]
                ys = [i0_v for _, i0_v in pts_list]
                self._active_marker_artist.set_data(xs, ys)
                self._active_marker_artist.set_visible(True)
            self.canvas_map.draw_idle()
        except Exception:
            if not self._map_redraw_timer.isActive():
                self._map_redraw_timer.start()

    def _on_scan_progress(self, prog: Dict[str, Any]) -> None:
        """Handle live intra-point turn progress and compute-node startup notification during a 2D scan."""
        if prog.get("event") == "queue_status":
            b_ids = ",".join(prog.get("batch_ids", []))
            st = prog.get("state", "PENDING")
            rsn = prog.get("reason", "Waiting")
            est = prog.get("start_est", "N/A")
            w_s = int(prog.get("wait_sec", 0))
            self.conn_badge.setText(f"Topaze Queue [{b_ids}]: {st} (Reason: {rsn}, Wait: {w_s}s)")
            self.conn_badge.setStyleSheet(
                "background: #fef3c7; color: #92400e; padding: 4px 8px; border-radius: 4px; font-weight: bold;"
            )
            self._set_status(
                f"[TGCC Queue {st} — {w_s}s] BatchID {b_ids} waiting for A100 slot | "
                f"SLURM Reason: {rsn} | Est. Start: {est} (Tip: switch to 'daemon_stream' for 0s wait)"
            )
            return

        if prog.get("event") == "job_started":
            node = prog.get("hostname", "a100_node")
            gpus = prog.get("gpus", [])
            s_i = int(prog.get("shard_index", 0))
            n_s = int(prog.get("num_shards", 1))
            w_k = int(prog.get("workers_per_gpu", self.f_workers_per_gpu.value()))
            g_str = f"{len(gpus)}x {gpus[0].get('name', 'A100')} [{gpus[0].get('backend', 'cuda')}]" if gpus else "A100"
            self.conn_badge.setText(f"Active Node: {node} [{g_str} × {w_k} Pts/GPU = {n_s * w_k} Parallel]")
            self.conn_badge.setStyleSheet(
                "background: #dcfce7; color: #166534; padding: 4px 8px; border-radius: 4px; font-weight: bold;"
            )
            self.set_status_msg = f"Running on A100 node {node} ({g_str}, {n_s * w_k} concurrent points)"
            self._set_status(self.set_status_msg)
            return

        now = time.perf_counter()
        if (now - self._last_progress_ui_ts) < 0.08:
            return
        self._last_progress_ui_ts = now

        pt_idx = int(prog.get("point_index", 1))
        tot_pts = max(1, int(prog.get("total_points", max(1, self._total_scan_points))))
        turn = int(prog.get("turn", 0))
        tot_turns = max(1, int(prog.get("total_turns", 1)))
        i0_mA = float(prog.get("I0_mA", 0.0))
        df_kHz = float(prog.get("Df_hc_kHz", 0.0))
        psi_deg = float(prog.get("psi_hc_deg", 0.0))
        xi_eff = float(prog.get("xi_eff", 0.0))
        dev = str(prog.get("device_name", "Server"))

        frac_overall = (self._completed_count + (turn / float(tot_turns))) / float(tot_pts)
        pct = min(99, max(0, int(round(100.0 * frac_overall))))
        self.progress_bar.setValue(pct)
        self._set_status(
            f"[2D Scan Pt {pt_idx}/{tot_pts} | {dev}] "
            f"I0={i0_mA:.0f}mA, Df={df_kHz:.1f}kHz (psi={psi_deg:.1f}°, xi={xi_eff:.2f}) — "
            f"Turn {turn}/{tot_turns} ({int(round(100*turn/tot_turns))}%)"
        )

        # Fast-update active marker positions on the 2D map without rebuilding figure
        new_marker = (df_kHz, i0_mA)
        prev_marker = self._active_scan_markers.get(pt_idx)
        self._active_scan_markers[pt_idx] = new_marker
        if prev_marker != new_marker:
            self._update_active_scan_marker_fast(new_marker)

        # Update live cavity waveform plot only if Tab 0 is currently visible
        if "hc_amp_kV" in prog and self.right_tabs.currentIndex() == 0:
            self._update_four_panel_waveforms(
                turn,
                prog.get("hc_amp_kV", []),
                prog.get("hc_phase_deg", []),
                prog.get("mc_amp_kV", []),
                prog.get("mc_phase_deg", []),
            )

    def _on_scan_point(self, pt: ScanPointTelemetry) -> None:
        """Immediately record completed point into 2D matrices and schedule coalesced Stability Map redraw."""
        self._completed_count += 1
        self._active_scan_markers.pop(int(pt.point_index), None)
        for k_id, (m_df, m_i0) in list(self._active_scan_markers.items()):
            if abs(m_df - pt.Df_hc_kHz) < 1e-3 and abs(m_i0 - pt.I0_mA) < 1e-3:
                self._active_scan_markers.pop(k_id, None)
        tot = max(1, pt.total_points)
        pct = int(round(100.0 * self._completed_count / tot))
        self.progress_bar.setValue(pct)

        if self._i0_grid_mA.size > 0 and self._df_grid_kHz.size > 0:
            i_idx = int(np.argmin(np.abs(self._i0_grid_mA - pt.I0_mA)))
            j_idx = int(np.argmin(np.abs(self._df_grid_kHz - pt.Df_hc_kHz)))
            if 0 <= i_idx < self._stability_code_mat.shape[0] and 0 <= j_idx < self._stability_code_mat.shape[1]:
                self._stability_code_mat[i_idx, j_idx] = float(pt.stability_code)
                self._growth_rate_mat[i_idx, j_idx] = float(pt.growth_rate_s1)
                if self._bunch_len_mat.shape == self._stability_code_mat.shape:
                    self._bunch_len_mat[i_idx, j_idx] = float(pt.mean_bunch_length_ps)
                    self._jitter_mat[i_idx, j_idx] = float(pt.centroid_jitter_ps)

        self._set_status(
            f"[2D Scan {self._completed_count}/{tot} ({pct}%) | GPU:{pt.gpu_id}] "
            f"I0={pt.I0_mA:.0f}mA, Df={pt.Df_hc_kHz:.1f}kHz (psi={pt.psi_hc_deg:.1f}°, xi={pt.xi_eff:.2f}) -> "
            f"{pt.status_label} (gr={pt.growth_rate_s1:+.1f} s^-1, {pt.elapsed_sec:.1f}s)"
        )
        self._append_log(
            f"[Point {self._completed_count}/{tot}] I0={pt.I0_mA:.0f}mA, Df={pt.Df_hc_kHz:.1f}kHz "
            f"(psi={pt.psi_hc_deg:.1f}°, xi={pt.xi_eff:.3f}) => {pt.status_label} (code={pt.stability_code}, gr={pt.growth_rate_s1:+.1f} s^-1)"
        )
        if not self._map_redraw_timer.isActive():
            self._map_redraw_timer.start()

    def _on_scan_completed(self, total_completed: int) -> None:
        self._set_top_busy(False)
        self._active_scan_marker = None
        self._active_scan_markers.clear()
        self.progress_bar.setValue(100)
        self._set_status(f"2D Stability Map Scan Done! ({total_completed} points completed)")
        self._append_log(f"[2D Scan Complete] Finished all {total_completed} points.")
        self._redraw_stability_map()
        self._prompt_save_settings_after_tracking("2D Stability Scan")

    # =========================================================================
    # Exact MATLAB `plotStabilityMap` Port (`STABLE_App.m` lines 4690-5038)
    # =========================================================================
    def _redraw_stability_map(self) -> None:
        self.fig_map.clf()
        self.ax_map = self.fig_map.add_subplot(1, 1, 1)

        use_psi = self.map_xaxis_combo.currentIndex() == 1
        cfg = self._collect_machine_config()

        if self._df_grid_kHz.size == 0 or self._i0_grid_mA.size == 0:
            s_cfg = self._collect_scan_config()
            Df_scan_kHz = np.array(s_cfg.get_df_grid_kHz(), dtype=np.float64)
            I0_scan_mA = np.array(s_cfg.get_i0_grid_mA(), dtype=np.float64)
            stability_code = np.full((len(I0_scan_mA), len(Df_scan_kHz)), np.nan)
            growth_rate = np.full((len(I0_scan_mA), len(Df_scan_kHz)), np.nan)
        else:
            Df_scan_kHz = self._df_grid_kHz
            I0_scan_mA = self._i0_grid_mA
            stability_code = self._stability_code_mat
            growth_rate = self._growth_rate_mat

        cspeed = 299792458.0
        f0_val = cspeed / cfg.C
        f_hc = cfg.n_hc * cfg.h * f0_val
        psi_scan_deg = np.rad2deg(np.arctan(2.0 * cfg.Q_hc * (Df_scan_kHz * 1e3) / f_hc))
        x_vals = psi_scan_deg if use_psi else Df_scan_kHz

        # 1. Background Growth-Rate Heatmap (STABLE_App.m lines 4712-4846)
        gr_plot = np.array(growth_rate, dtype=np.float64, copy=True)
        # Sanitize legacy placeholder values (e.g. 1/T0 ~ 8.46e5 from older scans, STABLE_App.m lines 4713-4719)
        legacy_mask = ~np.isnan(gr_plot) & (gr_plot > 5000.0)
        if np.any(legacy_mask):
            valid_pos = gr_plot[~legacy_mask & ~np.isnan(gr_plot) & (gr_plot > 0.0)]
            cap_val = float(max(200.0, 1.5 * np.max(valid_pos))) if valid_pos.size > 0 else 300.0
            gr_plot[legacy_mask] = cap_val

        gr_plot = np.where(np.isnan(gr_plot), np.nan, np.clip(gr_plot, -12.0, 600.0))

        all_complete = bool(gr_plot.shape[0] >= 2 and gr_plot.shape[1] >= 2 and np.all(~np.isnan(gr_plot)))
        if gr_plot.shape[0] >= 3 and gr_plot.shape[1] >= 3 and all_complete:
            k_sm = np.array([[0.05, 0.10, 0.05], [0.10, 0.40, 0.10], [0.05, 0.10, 0.05]])
            gr_plot = convolve(gr_plot, k_sm, mode="nearest")

        valid_gr = gr_plot[~np.isnan(gr_plot)]
        if valid_gr.size > 0 and np.max(valid_gr) > np.min(valid_gr):
            vmin, vmax = float(np.min(valid_gr)), float(max(np.max(valid_gr), 40.0))
        else:
            vmin, vmax = -40.0, 400.0

        # Smooth 400x300 cubic/makima interpolation when scan matrix is complete (STABLE_App.m lines 4775-4829)
        if all_complete and len(x_vals) >= 2 and len(I0_scan_mA) >= 2 and np.all(np.diff(x_vals) > 0) and np.all(np.diff(I0_scan_mA) > 0):
            kx = min(3, len(x_vals) - 1)
            ky = min(3, len(I0_scan_mA) - 1)
            spline = RectBivariateSpline(I0_scan_mA, x_vals, gr_plot, kx=ky, ky=kx)
            if use_psi:
                x_fine = np.linspace(float(np.min(x_vals)), float(np.max(x_vals)), 400)
                y_fine = np.linspace(float(np.min(I0_scan_mA)), float(np.max(I0_scan_mA)), 300)
                gr_fine = spline(y_fine, x_fine)
            else:
                d_df_f = float(Df_scan_kHz[1] - Df_scan_kHz[0]) / 2.0
                d_i0_f = float(I0_scan_mA[1] - I0_scan_mA[0]) / 2.0
                x_fine = np.linspace(float(np.min(Df_scan_kHz)) - d_df_f, float(np.max(Df_scan_kHz)) + d_df_f, 400)
                y_fine = np.linspace(float(np.min(I0_scan_mA)) - d_i0_f, float(np.max(I0_scan_mA)) + d_i0_f, 300)
                xq = np.clip(x_fine, float(np.min(Df_scan_kHz)), float(np.max(Df_scan_kHz)))
                yq = np.clip(y_fine, float(np.min(I0_scan_mA)), float(np.max(I0_scan_mA)))
                gr_fine = spline(yq, xq)
            mesh = self.ax_map.pcolormesh(
                x_fine, y_fine, gr_fine, cmap=_MATLAB_YLORRD_CMAP, vmin=vmin, vmax=vmax, shading="gouraud", zorder=1
            )
        else:
            if len(x_vals) > 1:
                dx = np.diff(x_vals)
                x_edges = np.concatenate(([x_vals[0] - dx[0] / 2.0], 0.5 * (x_vals[:-1] + x_vals[1:]), [x_vals[-1] + dx[-1] / 2.0]))
            else:
                x_edges = np.array([x_vals[0] - 0.5, x_vals[0] + 0.5])

            if len(I0_scan_mA) > 1:
                dy = np.diff(I0_scan_mA)
                y_edges = np.concatenate(([I0_scan_mA[0] - dy[0] / 2.0], 0.5 * (I0_scan_mA[:-1] + I0_scan_mA[1:]), [I0_scan_mA[-1] + dy[-1] / 2.0]))
            else:
                y_edges = np.array([I0_scan_mA[0] - 10.0, I0_scan_mA[0] + 10.0])

            masked_gr = np.ma.masked_invalid(gr_plot)
            mesh = self.ax_map.pcolormesh(
                x_edges, y_edges, masked_gr, cmap=_MATLAB_YLORRD_CMAP, vmin=vmin, vmax=vmax, shading="flat", zorder=1
            )
        cb = self.fig_map.colorbar(mesh, ax=self.ax_map)
        cb.set_label("Growth Rate  [$s^{-1}$]", fontsize=10)

        # 2. Analytical Reference Curves: Flat-potential xi_eff=1.0, xi_eff=1.05 & D-mode Thresh (STABLE_App.m lines 4847-4919)
        try:
            i0_min_line = float(np.min(I0_scan_mA)) if len(I0_scan_mA) > 1 else max(50.0, float(I0_scan_mA[0]) - 25.0)
            i0_max_line = float(np.max(I0_scan_mA)) if len(I0_scan_mA) > 1 else float(I0_scan_mA[0]) + 25.0
            cache_key = (
                round(i0_min_line, 2),
                round(i0_max_line, 2),
                round(cfg.V_mc, 1),
                round(cfg.U0, 1),
                round(cfg.E0, 1),
                round(cfg.C, 3),
                int(cfg.h),
                round(cfg.alpha_c, 8),
                round(cfg.sigma_e0, 7),
                int(cfg.n_hc),
                round(cfg.R_hc, 1),
                round(cfg.Q_hc, 1),
                bool(cfg.use_form_factor_iter),
            )
            if cache_key in self._ref_curves_cache:
                i0_line_A, df_100, df_105, psi_100, psi_105 = self._ref_curves_cache[cache_key]
            else:
                i0_line_A = np.linspace(i0_min_line, i0_max_line, 35) * 1e-3
                df_100, df_105, psi_100, psi_105 = [], [], [], []
                for ib in i0_line_A:
                    fp100 = calc_bunch_form_factor_eq(
                        "detune_for_xi", float(ib), cfg.V_mc, cfg.U0, cfg.E0, cfg.C, cfg.h,
                        cfg.alpha_c, cfg.sigma_e0, cfg.n_hc, cfg.R_hc, cfg.Q_hc, 1.0, cfg.use_form_factor_iter,
                    )
                    fp105 = calc_bunch_form_factor_eq(
                        "detune_for_xi", float(ib), cfg.V_mc, cfg.U0, cfg.E0, cfg.C, cfg.h,
                        cfg.alpha_c, cfg.sigma_e0, cfg.n_hc, cfg.R_hc, cfg.Q_hc, 1.05, cfg.use_form_factor_iter,
                    )
                    df_100.append(fp100["fre_shift_hc"] * 1e-3)
                    df_105.append(fp105["fre_shift_hc"] * 1e-3)
                    psi_100.append(fp100["psi_hc_deg"])
                    psi_105.append(fp105["psi_hc_deg"])
                self._ref_curves_cache[cache_key] = (i0_line_A, df_100, df_105, psi_100, psi_105)

            lbl_100 = r"Flat-pot. ($\xi_{\mathrm{eff}}=1.0, F_n$)" if cfg.use_form_factor_iter else r"Flat-pot. ($\xi=1.0$)"
            lbl_105 = r"$\xi_{\mathrm{eff}}=1.05$ (PTBL, $F_n$)" if cfg.use_form_factor_iter else r"$\xi=1.05$ (PTBL)"

            if use_psi:
                self.ax_map.plot(psi_100, i0_line_A * 1e3, "k--", lw=2.0, label=lbl_100, zorder=3)
                self.ax_map.plot(psi_105, i0_line_A * 1e3, "g:", lw=1.8, label=lbl_105, zorder=3)
            else:
                self.ax_map.plot(df_100, i0_line_A * 1e3, "k--", lw=2.0, label=lbl_100, zorder=3)
                self.ax_map.plot(df_105, i0_line_A * 1e3, "g:", lw=1.8, label=lbl_105, zorder=3)
        except Exception:
            pass

        # 3. Fitted Physical Safety Boundary (STABLE_App.m lines 4921-4956 & fitAndPlotSafetyBoundary lines 2030-2052)
        try:
            fit_res = fit_safety_boundary(Df_scan_kHz, I0_scan_mA, stability_code)
            if fit_res is not None:
                i0_dense = fit_res["i0_dense"]
                df_dense = fit_res["df_dense"]
                r2_val = fit_res["r2"]
                if use_psi:
                    psi_dense = np.rad2deg(np.arctan(2.0 * cfg.Q_hc * (df_dense * 1e3) / f_hc))
                    lbl_sb = f"Safety Boundary ($\\psi_{{\\mathrm{{crit}}}}$, $R^2={r2_val:.3f}$)" if self._show_critical_boundary_dots else r"Safety Boundary ($\psi_{\mathrm{crit}}$)"
                    self.ax_map.plot(psi_dense, i0_dense, "-", color=(0.05, 0.30, 0.85), lw=2.5, label=lbl_sb, zorder=3)
                    if self._show_critical_boundary_dots:
                        psi_crit_pts = np.rad2deg(np.arctan(2.0 * cfg.Q_hc * (fit_res["crit_df"] * 1e3) / f_hc))
                        self.ax_map.plot(psi_crit_pts, fit_res["crit_i0"], "k.", markersize=12, label="Critical Trans. Points", zorder=6)
                else:
                    a_v, b_v, sgn = fit_res["a"], abs(fit_res["b"]), fit_res["sgn"]
                    lbl_sb = (
                        f"Safety: $\\Delta f > {a_v:.3f} I_0 {sgn} {b_v:.1f}$ ($R^2={r2_val:.3f}$)"
                        if self._show_critical_boundary_dots
                        else f"Safety: $\\Delta f > {a_v:.3f} I_0 {sgn} {b_v:.1f}$"
                    )
                    self.ax_map.plot(df_dense, i0_dense, "-", color=(0.05, 0.30, 0.85), lw=2.5, label=lbl_sb, zorder=3)
                    if self._show_critical_boundary_dots:
                        self.ax_map.plot(fit_res["crit_df"], fit_res["crit_i0"], "k.", markersize=12, label="Critical Trans. Points", zorder=6)
        except Exception:
            pass

        # 4. Simulation Point Markers matching MATLAB STABLE_App.m (lines 4958-4983)
        X_grid, I0_grid = np.meshgrid(x_vals, I0_scan_mA)
        m0 = stability_code == 0
        m1 = stability_code == 1
        m2 = stability_code == 2
        m3 = stability_code == 3
        m4 = stability_code == 4

        # Proxy legend entries (always visible like STABLE_App.m lines 4959-4963)
        self.ax_map.plot([], [], "o", markerfacecolor=(0.18, 0.72, 0.22), markeredgecolor=(0.0, 0.28, 0.0), markersize=7, label="Stable (o)")
        self.ax_map.plot([], [], "x", color="k", markeredgewidth=2.2, markersize=8, label="PTBL (x)")
        self.ax_map.plot([], [], "*", markerfacecolor=(0.10, 0.55, 0.98), markeredgecolor="k", markersize=10, label="Coupled D-Q (\u2605)")
        self.ax_map.plot([], [], "^", markerfacecolor=(0.58, 0.12, 0.78), markeredgecolor="k", markersize=7, label="LCBI (\u25b2)")
        self.ax_map.plot([], [], "D", markerfacecolor=(0.45, 0.45, 0.45), markeredgecolor="k", markersize=6, label="Other (\u25c6)")

        if np.any(m0):
            self.ax_map.plot(X_grid[m0], I0_grid[m0], "o", markerfacecolor=(0.18, 0.72, 0.22), markeredgecolor=(0.0, 0.28, 0.0), markersize=6.5, zorder=4)
        if np.any(m1):
            self.ax_map.plot(X_grid[m1], I0_grid[m1], "x", color="w", markeredgewidth=3.6, markersize=9, zorder=4)
            self.ax_map.plot(X_grid[m1], I0_grid[m1], "x", color="k", markeredgewidth=2.0, markersize=8, zorder=5)
        if np.any(m2):
            self.ax_map.plot(X_grid[m2], I0_grid[m2], "*", markerfacecolor=(0.10, 0.55, 0.98), markeredgecolor="k", markersize=10, zorder=4)
        if np.any(m3):
            self.ax_map.plot(X_grid[m3], I0_grid[m3], "^", markerfacecolor=(0.58, 0.12, 0.78), markeredgecolor="k", markersize=7.5, zorder=4)
        if np.any(m4):
            self.ax_map.plot(X_grid[m4], I0_grid[m4], "D", markerfacecolor=(0.50, 0.50, 0.50), markeredgecolor="k", markersize=6.5, zorder=4)

        # Active scanning point indicator(s) (red hollow rings for all concurrent points in flight)
        pts_list = list(self._active_scan_markers.values())
        if not pts_list and self._active_scan_marker is not None:
            pts_list = [self._active_scan_marker]
        if pts_list:
            xs = [
                (math.degrees(math.atan(2.0 * cfg.Q_hc * (df_v * 1e3) / f_hc)) if use_psi else df_v)
                for df_v, _ in pts_list
            ]
            ys = [i0_v for _, i0_v in pts_list]
            (self._active_marker_artist,) = self.ax_map.plot(
                xs, ys, "o", markerfacecolor="none", markeredgecolor="#ef4444",
                markeredgewidth=2.5, markersize=12, label="Tracking Now", zorder=6
            )
        else:
            (self._active_marker_artist,) = self.ax_map.plot(
                [], [], "o", markerfacecolor="none", markeredgecolor="#ef4444",
                markeredgewidth=2.5, markersize=12, zorder=6
            )

        # 5. Titles, Parameter Summary Subtitle & Axis Limits (STABLE_App.m lines 4985-5030)
        pi_str = f"PI(K_p={cfg.pi_kp:.2g}, K_i={cfg.pi_ki:.0e})" if cfg.pi_enable else "PI: Off"
        if not cfg.dfb_enable:
            dfb_str = "DFB: Off"
        elif cfg.dfb_adaptive:
            dfb_str = "DFB: Adaptive"
        elif cfg.dfb_auto_scale:
            dfb_str = f"DFB: Auto({cfg.dfb_G1:.2g}-{cfg.dfb_G2:.2g})"
        else:
            dfb_str = f"DFB(g={cfg.dfb_gain:.2g})"
        hom_str = ("HOM: Aligned" if cfg.hom_auto_align else "HOM: Nominal") if cfg.homs else "HOM: Off"
        param_sub = (
            f"V_mc={cfg.V_mc*1e-6:.2f} MV, R_hc={cfg.R_hc*1e-6:.2f} MOhm, Q_hc={cfg.Q_hc:.0f}  |  "
            f"{pi_str}  |  {dfb_str}  |  {hom_str}"
        )

        if use_psi:
            self.ax_map.set_title(f"2D Stability Map vs HC Tuning Angle $\\psi_{{HC}}$\n{param_sub}", fontsize=10)
            self.ax_map.set_xlabel("HHC Tuning Angle $\\psi_{HC}$  [$^\\circ$]", fontsize=10.5)
            p_min, p_max = float(np.min(psi_scan_deg)), float(np.max(psi_scan_deg))
            if (p_max - p_min) < 1.5:
                # High-Q superconducting cavity (e.g. Q_hc ~ 2e8 where psi_HC is near 89.99 deg)
                pad_p = max(1e-4, 0.08 * (p_max - p_min))
                self.ax_map.set_xlim(p_min - pad_p, p_max + pad_p)
            else:
                psi_lo = min(65.0, math.floor(p_min) - 0.5)
                psi_hi = min(90.0, max(85.0, math.ceil(p_max * 2.0) / 2.0 + 0.2))
                self.ax_map.set_xlim(psi_lo, psi_hi)
        else:
            self.ax_map.set_title(f"2D Stability Map vs $\\Delta f_{{HC}}$\n{param_sub}", fontsize=10)
            self.ax_map.set_xlabel("HHC Detuning $\\Delta f_{HC}$  [kHz]", fontsize=10.5)
            d_df = float(Df_scan_kHz[1] - Df_scan_kHz[0]) / 2.0 if len(Df_scan_kHz) > 1 else 2.5
            self.ax_map.set_xlim(float(np.min(Df_scan_kHz)) - d_df, float(np.max(Df_scan_kHz)) + d_df)

        self.ax_map.set_ylabel("Beam Current $I_0$  [mA]", fontsize=10.5)
        d_i0 = float(I0_scan_mA[1] - I0_scan_mA[0]) / 2.0 if len(I0_scan_mA) > 1 else 12.5
        self.ax_map.set_ylim(float(np.min(I0_scan_mA)) - d_i0, float(np.max(I0_scan_mA)) + d_i0)

        self.ax_map.grid(True, linestyle=":", alpha=0.55)
        loc_map = {
            "Lower Right": "lower right",
            "Upper Right": "upper right",
            "Upper Left": "upper left",
            "Lower Left": "lower left",
            "Best (Auto)": "best",
        }
        leg_choice = self.map_legend_combo.currentText() if hasattr(self, "map_legend_combo") else "Lower Right"
        if leg_choice != "Hide":
            leg = self.ax_map.legend(loc=loc_map.get(leg_choice, "lower right"), fontsize=8.5, framealpha=0.88)
            if leg is not None:
                leg.set_draggable(True)
        self.canvas_map.draw_idle()

    # =========================================================================
    # Offline Scan Analysis Tab (`AnalyzeButtonPushed` & `plotAnalysisSelection` Port)
    # =========================================================================
    def _on_analyze_scan_folder(self) -> None:
        default_dir = str(Path(__file__).resolve().parents[3] / "scan_data")
        folder = QFileDialog.getExistingDirectory(
            self, "Select Scan Data Folder with SOLEIL_I0*detune*.mat", default_dir
        )
        if not folder:
            return
        try:
            sum_f = Path(folder) / "scan_summary.mat"
            if sum_f.exists():
                self._load_scan_summary_path(sum_f, restore_config=True, switch_tab=False)
            self._offline_analysis_data = analyze_scan_folder_mat(folder)
            self.right_tabs.setCurrentIndex(2)
            self._redraw_offline_analysis_plot()
            best_df = self._offline_analysis_data["best_detune_kHz"]
            best_bl = self._offline_analysis_data["best_bunch_length_ps"]
            if not math.isnan(best_df):
                self._set_status(f"Optimal: {best_df:.0f} kHz (BL: {best_bl:.2f} ps)")
        except Exception as exc:
            QMessageBox.warning(self, "Analysis Error", str(exc))

    def _on_plot_current_scan_slice(self) -> None:
        """Plot Bunch Length & Centroid Jitter vs Df_HC from the in-memory 2D scan."""
        if self._df_grid_kHz.size == 0 or self._bunch_len_mat.size == 0 or np.all(np.isnan(self._bunch_len_mat)):
            QMessageBox.information(
                self,
                "Offline Scan Analysis",
                "Run a 2D scan in this session first, or click 'Analyze Scan Folder (SOLEIL_I0*.mat)...' to analyze saved MATLAB point files.",
            )
            return
        # Pick row with most valid points
        valid_counts = np.sum(~np.isnan(self._bunch_len_mat), axis=1)
        r_idx = int(np.argmax(valid_counts))
        valid = ~np.isnan(self._bunch_len_mat[r_idx, :])
        x_df = self._df_grid_kHz[valid]
        y_bl = self._bunch_len_mat[r_idx, valid]
        y_jit = self._jitter_mat[r_idx, valid]
        y_es = np.full_like(y_bl, self.f_sigma_e0.value())
        stable_idx = np.where(y_jit < 1.5)[0]
        if stable_idx.size > 0:
            b_i = stable_idx[int(np.argmax(y_bl[stable_idx]))]
            best_df, best_bl = float(x_df[b_i]), float(y_bl[b_i])
        else:
            b_i = int(np.argmin(y_jit)) if y_jit.size > 0 else 0
            best_df = float(x_df[b_i]) if x_df.size > 0 else float("nan")
            best_bl = float(y_bl[b_i]) if y_bl.size > 0 else float("nan")
        self._offline_analysis_data = {
            "scan_kHz": x_df,
            "bunch_length_ps": y_bl,
            "centroid_mean_ps": y_jit,
            "centroid_jitter_ps": y_jit,
            "energy_spread": y_es,
            "best_detune_kHz": best_df,
            "best_bunch_length_ps": best_bl,
        }
        self._redraw_offline_analysis_plot()

    def _redraw_offline_analysis_plot(self) -> None:
        """Exact port of MATLAB `plotAnalysisSelection` (STABLE_App.m lines 1418-1557)."""
        self.fig_an.clf()
        self.ax_an = self.fig_an.add_subplot(1, 1, 1)
        sel = self.analysis_mode_combo.currentText()

        if sel == "Mode-0 Spectrum (D-mode Check)":
            fr = self._last_telemetry_frame
            if fr is None or not fr.turns_axis or len(fr.centroid_jitter_ps_history) < 4:
                self.ax_an.set_title("Run a single-point tracking first to inspect Mode-0 Coherent Spectrum")
                self.canvas_an.draw_idle()
                return
            cfg = self._collect_machine_config()
            dt_rec = cfg.Recor_step * (cfg.C / 299792458.0)
            sig = np.asarray(fr.centroid_jitter_ps_history, dtype=np.float64)
            sig = sig - np.mean(sig)
            freqs_kHz = np.fft.rfftfreq(len(sig), d=dt_rec) * 1e-3
            amp = np.abs(np.fft.rfft(sig))
            amp_norm = amp / max(1e-12, float(np.max(amp)))
            self.ax_an.plot(freqs_kHz, amp_norm, "-", lw=1.8, color="#1973c7", label=r"STABLE $\langle\tau\rangle_0$ FFT")
            self.ax_an.set_xlabel("Coherent Frequency $f$ [kHz]")
            self.ax_an.set_ylabel("Normalized Spectral Amplitude")
            self.ax_an.set_title(f"Mode-0 Coherent Spectrum (Dominant f = {fr.dominant_freq_Hz*1e-3:.2f} kHz)")
            self.ax_an.set_ylim(0.0, 1.12)
            self.ax_an.grid(True, alpha=0.35)
            self.ax_an.legend(loc="upper right")
            self.canvas_an.draw_idle()
            return

        if self._offline_analysis_data is None or self._offline_analysis_data["scan_kHz"].size == 0:
            self.ax_an.set_title("Click 'Analyze Scan Folder (SOLEIL_I0*.mat)...' or 'Plot Current 2D Scan Slice'")
            self.ax_an.set_xlabel(r"HHC Detuning $\Delta f_{HC}$ [kHz]")
            self.ax_an.grid(True, alpha=0.35)
            self.canvas_an.draw_idle()
            return

        d = self._offline_analysis_data
        x = d["scan_kHz"]
        y_bl = d["bunch_length_ps"]
        y_c = d["centroid_jitter_ps"]
        y_e = d["energy_spread"]
        best_df = d["best_detune_kHz"]
        best_bl = d["best_bunch_length_ps"]

        self.ax_an.set_title("Bunch Lengthening Effect vs. HHC Detuning")
        self.ax_an.set_xlabel(r"HHC Detuning $\Delta f_{HC}$ [kHz]")
        self.ax_an.grid(True, alpha=0.35)

        c_blue = (0.0, 0.4470, 0.7410)
        c_orange = (0.8500, 0.3250, 0.0980)
        c_green = (0.4660, 0.6740, 0.1880)

        if sel == "Bunch Length":
            self.ax_an.plot(x, y_bl, "-o", lw=2.0, color=c_blue, label="Bunch length")
            self.ax_an.set_ylabel(r"Bunch Length $\sigma_\tau$ [ps]", color=c_blue)
            self.ax_an.legend(loc="best")
        elif sel == "Centroid":
            self.ax_an.plot(x, y_c, "--s", lw=1.5, color=c_orange, label="Centroid jitter")
            self.ax_an.set_ylabel("Centroid [ps]", color=c_orange)
            self.ax_an.legend(loc="best")
        elif sel == "Energy Spread":
            self.ax_an.plot(x, y_e, "-.d", lw=1.5, color=c_green, label="Energy spread")
            self.ax_an.set_ylabel(r"Energy spread $\sigma_\delta$", color=c_green)
            self.ax_an.legend(loc="best")
        elif sel == "Bunch Length + Centroid":
            l1 = self.ax_an.plot(x, y_bl, "-o", lw=2.0, color=c_blue, label="Bunch length")
            self.ax_an.set_ylabel(r"Bunch Length $\sigma_\tau$ [ps]", color=c_blue)
            ax_r = self.ax_an.twinx()
            l2 = ax_r.plot(x, y_c, "--s", lw=1.5, color=c_orange, label="Centroid")
            ax_r.set_ylabel("Centroid [ps]", color=c_orange)
            lines = l1 + l2
            self.ax_an.legend(lines, [l.get_label() for l in lines], loc="best")
        elif sel == "Bunch Length + Energy Spread":
            l1 = self.ax_an.plot(x, y_bl, "-o", lw=2.0, color=c_blue, label="Bunch length")
            self.ax_an.set_ylabel(r"Bunch Length $\sigma_\tau$ [ps]", color=c_blue)
            ax_r = self.ax_an.twinx()
            l2 = ax_r.plot(x, y_e, "-.d", lw=1.5, color=c_green, label="Energy spread")
            ax_r.set_ylabel(r"Energy spread $\sigma_\delta$", color=c_green)
            lines = l1 + l2
            self.ax_an.legend(lines, [l.get_label() for l in lines], loc="best")
        elif sel == "Centroid + Energy Spread":
            l1 = self.ax_an.plot(x, y_c, "--s", lw=1.5, color=c_orange, label="Centroid")
            self.ax_an.set_ylabel("Centroid [ps]", color=c_orange)
            ax_r = self.ax_an.twinx()
            l2 = ax_r.plot(x, y_e, "-.d", lw=1.5, color=c_green, label="Energy spread")
            ax_r.set_ylabel(r"Energy spread $\sigma_\delta$", color=c_green)
            lines = l1 + l2
            self.ax_an.legend(lines, [l.get_label() for l in lines], loc="best")
        else:
            l1 = self.ax_an.plot(x, y_bl, "-o", lw=2.0, color=c_blue, label="Bunch length")
            self.ax_an.set_ylabel(r"Bunch Length $\sigma_\tau$ [ps]", color=c_blue)
            ax_r = self.ax_an.twinx()
            l2 = ax_r.plot(x, y_c, "--s", lw=1.5, color=c_orange, label="Centroid")
            l3 = ax_r.plot(x, y_e, "-.d", lw=1.5, color=c_green, label="Energy spread")
            ax_r.axhline(1.5, color="r", ls=":", lw=1.2, label="Stability Limit (1.5 ps)")
            ax_r.set_ylabel("Centroid [ps] / Energy spread", color=c_orange)
            lines = l1 + l2 + l3
            self.ax_an.legend(lines, [l.get_label() for l in lines], loc="best")

        if not math.isnan(best_df) and sel in ("Bunch Length", "Bunch Length + Centroid", "Bunch Length + Energy Spread", "All"):
            self.ax_an.plot([best_df], [best_bl], "*", ms=14, mfc="r", mec="r", zorder=6)
            dx_off = max(1.0, (float(np.max(x)) - float(np.min(x))) * 0.04) if x.size > 1 else 2.0
            self.ax_an.text(best_df + dx_off, best_bl, f"Optimal ({best_df:.0f} kHz)", color="r", fontsize=10, fontweight="bold")

        self.canvas_an.draw_idle()

    def _refresh_nyquist_plot(self) -> None:
        cfg = self._collect_machine_config()
        try:
            # Compute directly in memory (< 2 ms) matching MATLAB updateLLRFNyquistPlot (STABLE_App.m lines 1245-1408)
            cfg.sync_derived_cavity_params()
            df_mc = (
                compute_optimal_mc_detune(
                    cfg.C, cfg.h, cfg.I0, cfg.V_mc, cfg.U0, cfg.R_mc, cfg.Q_mc,
                    cfg.R_hc, cfg.Q_hc, cfg.fre_shift_hc, cfg.n_hc,
                )
                if cfg.mc_auto_detune
                else cfg.fre_shift_mc
            )
            inc_detune = self.chk_llrf_detune.isChecked() if hasattr(self, "chk_llrf_detune") else True
            loop_mode = self.llrf_mode_combo.currentText() if hasattr(self, "llrf_mode_combo") else "Combined (DFB + PI)"

            use_pi = cfg.pi_enable
            use_dfb = cfg.dfb_enable
            if loop_mode == "Direct RF Feedback (DFB Only)":
                use_pi, use_dfb = False, True
            elif loop_mode == "PI Feedback (PI Only)":
                use_pi, use_dfb = True, False
            elif loop_mode == "Combined (DFB + PI)":
                use_pi, use_dfb = True, True

            for ax in (self.ax_bode_mag, self.ax_bode_pha, self.ax_nyquist):
                ax.clear()

            nyq = analyze_llrf_nyquist(
                C=cfg.C,
                h=cfg.h,
                Q_mc=cfg.Q_mc,
                fre_shift_mc=df_mc,
                pi_enable=use_pi,
                pi_kp=cfg.pi_kp,
                pi_ki=cfg.pi_ki,
                pi_m=cfg.pi_m,
                pi_d=cfg.pi_d,
                dfb_enable=use_dfb,
                dfb_gain=cfg.get_dfb_gain_for_current(cfg.I0),
                dfb_phase_shift_deg=cfg.dfb_phase_shift_deg,
                dfb_delay=cfg.dfb_delay,
                dfb_m=cfg.dfb_m,
                include_detune=inc_detune,
            )
            f_kHz = np.array(nyq["f_pos_kHz"])

            if loop_mode == "Gain Sweep (Compare G_dfb)":
                g_max_val = float(nyq["g_dfb_max"])
                g_Sweep = [0.15, 0.25, 0.50, min(0.95, g_max_val if math.isfinite(g_max_val) else 0.85)]
                colors = ["#0284c7", "#16a34a", "#d97706", "#dc2626"]
                for g_v, col in zip(g_Sweep, colors):
                    n_g = analyze_llrf_nyquist(
                        C=cfg.C, h=cfg.h, Q_mc=cfg.Q_mc, fre_shift_mc=df_mc,
                        pi_enable=cfg.pi_enable, pi_kp=cfg.pi_kp, pi_ki=cfg.pi_ki, pi_m=cfg.pi_m, pi_d=cfg.pi_d,
                        dfb_enable=True, dfb_gain=g_v, dfb_phase_shift_deg=cfg.dfb_phase_shift_deg,
                        dfb_delay=cfg.dfb_delay, dfb_m=cfg.dfb_m, include_detune=inc_detune,
                    )
                    lbl_g = f"G_dfb={g_v:.2f} (GM={n_g['gain_margin_dB']:.1f}dB)"
                    self.ax_bode_mag.semilogx(f_kHz, n_g["mag_dB"], color=col, lw=1.5, label=lbl_g)
                    self.ax_bode_pha.semilogx(f_kHz, n_g["phase_deg"], color=col, lw=1.5)
                    self.ax_nyquist.plot(n_g["nyq_pos_real"], n_g["nyq_pos_imag"], color=col, lw=1.6, label=f"G={g_v:.2f}")
                self.ax_bode_mag.legend(loc="lower left", fontsize=7.5)
            else:
                mag_dB = np.array(nyq["mag_dB"])
                pha_deg = np.array(nyq["phase_deg"])
                re_p = np.array(nyq["nyq_pos_real"])
                im_p = np.array(nyq["nyq_pos_imag"])
                re_n = np.array(nyq["nyq_neg_real"])
                im_n = np.array(nyq["nyq_neg_imag"])
                self.ax_bode_mag.semilogx(f_kHz, mag_dB, color="#2563eb", lw=1.7, label=loop_mode)
                self.ax_bode_pha.semilogx(f_kHz, pha_deg, color="#2563eb", lw=1.7)
                self.ax_nyquist.plot(re_p, im_p, color="#2563eb", lw=1.8, label="+w_m branch")
                self.ax_nyquist.plot(re_n, im_n, color="#2563eb", ls=":", lw=1.1, label="-w_m branch")

            self.ax_bode_mag.axhline(0.0, color="gray", ls="--", lw=0.9)
            self.ax_bode_mag.set_ylabel("Magnitude [dB]")
            self.ax_bode_mag.set_ylim(-45.0, 25.0)
            self.ax_bode_mag.set_title(f"Open-Loop Bode (GM={nyq['gain_margin_dB']:.1f} dB, PM={nyq['phase_margin_deg']:.1f}°)")
            self.ax_bode_mag.grid(True, which="both", alpha=0.35)

            self.ax_bode_pha.axhline(-180.0, color="#dc2626", ls="--", lw=1.2)
            if math.isfinite(nyq.get("f_pi_kHz", float("nan"))):
                self.ax_bode_pha.plot([nyq["f_pi_kHz"]], [-180.0], "ro", ms=5)
            self.ax_bode_pha.set_xlabel("Modulation Frequency f_m [kHz]")
            self.ax_bode_pha.set_ylabel("Phase [deg]")
            self.ax_bode_pha.set_ylim(-360.0, 25.0)
            self.ax_bode_pha.grid(True, which="both", alpha=0.35)

            th = np.linspace(0.0, 2.0 * math.pi, 240)
            self.ax_nyquist.plot(np.cos(th), np.sin(th), color="gray", ls="--", lw=1.0, label="Unit Circle")
            self.ax_nyquist.axvline(0.0, color="gray", ls=":", lw=0.8)
            self.ax_nyquist.axhline(0.0, color="gray", ls=":", lw=0.8)
            self.ax_nyquist.plot(-1.0, 0.0, "r+", ms=12, mew=2.2, label="Critical (-1, j0)")
            self.ax_nyquist.plot(-1.0, 0.0, "ro", ms=7, mfc="none", mew=1.4)
            self.ax_nyquist.set_xlim(-2.2, 2.5)
            self.ax_nyquist.set_ylim(-2.2, 2.2)
            self.ax_nyquist.set_aspect("equal", adjustable="box")
            self.ax_nyquist.set_xlabel("Re[L(jw_m)]")
            self.ax_nyquist.set_ylabel("Im[L(jw_m)]")
            self.ax_nyquist.set_title(f"LLRF Nyquist (Max G_dfb = {nyq['g_dfb_max']:.2f})")
            self.ax_nyquist.legend(loc="upper right", fontsize=8)
            self.ax_nyquist.grid(True, alpha=0.35)

            if hasattr(self, "lbl_llrf_metrics"):
                gm_v = float(nyq["gain_margin_dB"])
                pm_v = float(nyq["phase_margin_deg"])
                f_pi = float(nyq.get("f_pi_kHz", float("nan")))
                g_max = float(nyq["g_dfb_max"])
                is_stable = gm_v > 0.0 and (math.isnan(pm_v) or pm_v > 0.0)
                st_word = "STABLE" if is_stable else "UNSTABLE (Encircles -1+j0)"
                pm_s = f"{pm_v:.1f} deg" if math.isfinite(pm_v) else "Inf"
                self.lbl_llrf_metrics.setText(
                    f"GM = {gm_v:.1f} dB (@f_pi={f_pi:.1f} kHz)  |  PM = {pm_s}  |  Max G_dfb = {g_max:.2f}  |  {st_word}"
                )
                self.lbl_llrf_metrics.setStyleSheet(
                    f"color: {'#15803d' if is_stable else '#dc2626'}; font-weight: bold;"
                )

            self.canvas_nyq.draw_idle()
        except Exception as exc:
            self._append_log(f"[Nyquist Error] {exc}")

    def _on_job_error(self, err: str) -> None:
        self._set_top_busy(False)
        self._set_status(f"Error: {err}")
        self._append_log(f"[ERROR] {err}")

    def _append_log(self, msg: str) -> None:
        self.log_console.append(msg)

    def _set_status(self, txt: str) -> None:
        self.status_label.setText(txt)

    def closeEvent(self, event) -> None:  # noqa: N802
        try:
            self.vault.save_last_session_state(
                self.session,
                self._collect_machine_config(),
                self._collect_scan_config(),
            )
        except Exception:
            pass
        self.connector.disconnect()
        super().closeEvent(event)
