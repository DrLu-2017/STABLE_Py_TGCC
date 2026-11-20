"""Core types, data structures, and physical constants for the STABLE acceleration physics framework."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import IntEnum
from typing import NamedTuple, Optional
import numpy as np


class StabilityStatus(IntEnum):
    """Longitudinal multi-bunch instability status classification."""
    STABLE = 0       # Completely stable (Landau-damped flat-potential island or single-well damped)
    PTBL = 1         # Periodic Transient Bunch Lengthening / Sawtooth relaxation instability
    COUPLED_DQ = 2   # Coupled Dipole-Quadrupole mode oscillation (transition regime)
    LCBI = 3         # Longitudinal Coupled-Bunch Instability (HOM-driven)
    LOST = 4         # Beam loss, extreme blowup, or numerical divergence


@dataclass(slots=True, frozen=True)
class InstabilityResult:
    """Standardized output of the instability evaluation pipeline."""
    status: StabilityStatus
    is_stable: bool
    mode_id: int                    # Dominant coupled-bunch mode mu in [0, h-1]
    growth_rate: float              # Exponential growth rate [s^-1]
    bunch_length_factor: float      # Relative bunch lengthening factor: sigma_z / sigma_z0
    centroid_jitter_ps: float       # RMS bunch centroid jitter [ps]
    energy_spread_factor: float     # Relative energy spread: sigma_delta / sigma_delta0
    crest_factor: float             # Spatial DFT mode power crest factor
    nyquist_margin: float = 0.0     # LLRF Nyquist stability / distance margin
    details: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "status": self.status.name,
            "status_code": int(self.status),
            "is_stable": self.is_stable,
            "mode_id": self.mode_id,
            "growth_rate": self.growth_rate,
            "bunch_length_factor": self.bunch_length_factor,
            "centroid_jitter_ps": self.centroid_jitter_ps,
            "energy_spread_factor": self.energy_spread_factor,
            "crest_factor": self.crest_factor,
            "nyquist_margin": self.nyquist_margin,
        }


@dataclass(slots=True)
class MachineParameters:
    """Storage ring accelerator and RF system configuration parameters (SOLEIL II IDs-open worst-case RF default)."""
    # Ring geometry and beam energy
    C: float = 354.097                  # Circumference [m] (f_rf = 352.20 MHz)
    E0: float = 2.75e9                  # Nominal beam energy [eV]
    U0: float = 470e3                   # Energy loss per turn [eV] (IDs open)
    h: int = 416                        # Harmonic number
    alpha_c: float = 1.07e-4            # Momentum compaction factor
    tau_z: float = 13.7e-3              # Longitudinal radiation damping time [s] (IDs open, weakest damping)
    sigma_delta0: float = 0.92e-3       # Natural RMS energy spread (IDs open, weakest Landau damping)
    sigma_t0: float = 9.02e-12          # Natural RMS bunch length [s]

    # Main Cavity (MC, 4 cavities, beta = 4.95 -> Q_L = 6000, R_L = 3.36 MOhm)
    V_mc: float = 1.7e6                 # Main cavity total peak voltage [V]
    Q_mc: float = 6000.0                # Loaded quality factor of MC
    R_mc: float = 3.36e6                # Total loaded shunt impedance of MC [Ohm]
    RoQ_mc: float = 560.0               # Total geometric R/Q of 4 MCs [Ohm] (4 * 140 Ohm)
    fre_shift_mc: float = -35e3         # MC detuning frequency [Hz]

    # Harmonic Cavity (HC, 2 passive normal-conducting copper cavities)
    n_hc: int = 4                       # Harmonic ratio (4th harmonic for SOLEIL-II)
    Q_hc: float = 31000.0               # Loaded quality factor of HC
    R_hc: float = 1.86e6                # Total loaded shunt impedance of HC [Ohm]
    RoQ_hc: float = 60.0                # Total geometric R/Q of 2 HCs [Ohm] (2 * 30 Ohm)
    fre_shift_hc: float = 70e3          # HC detuning frequency [Hz]
    use_form_factor: bool = True        # Toggle iterative stretched bunch form factor F_n(sigma_tau)

    # Beam current and filling
    I0: float = 0.338                   # Total average beam current [A]

    # Direct RF Feedback (DFB) & LLRF
    DFB_enable: bool = True
    DFB_gain: float = 0.20              # DFB proportional gain
    DFB_phase_shift: float = 0.0        # DFB phase shift [rad]
    DFB_delay_buckets: int = 704        # Feedback loop delay in RF buckets (2.0 us = 704 * Tb)

    # Physical constants
    c_speed: float = 299792458.0

    @property
    def T0(self) -> float:
        """Revolution period [s]."""
        return self.C / self.c_speed

    @property
    def f0(self) -> float:
        """Revolution frequency [Hz]."""
        return 1.0 / self.T0

    @property
    def omega0(self) -> float:
        """Revolution angular frequency [rad/s]."""
        return 2.0 * np.pi * self.f0

    @property
    def frf(self) -> float:
        """Fundamental RF frequency [Hz]."""
        return self.h * self.f0

    @property
    def fhc(self) -> float:
        """Harmonic cavity resonant/reference frequency [Hz]."""
        return self.n_hc * self.frf

    @property
    def damping_rate(self) -> float:
        """Radiation damping rate 1 / tau_z [s^-1]."""
        return 1.0 / self.tau_z

    @property
    def c(self) -> float:
        """Speed of light [m/s]."""
        return self.c_speed

    @property
    def f_rf(self) -> float:
        """Alias for fundamental RF frequency [Hz]."""
        return self.frf

    @property
    def sigma_delta(self) -> float:
        """Alias for natural RMS energy spread."""
        return self.sigma_delta0

    @property
    def R_over_Q_hc(self) -> float:
        """Alias for R/Q of harmonic cavity [Ohm]."""
        return self.RoQ_hc

    @property
    def dfb_gain(self) -> float:
        """Alias for DFB gain."""
        return self.DFB_gain

    @property
    def dfb_delay(self) -> float:
        """Alias for DFB delay in seconds."""
        return self.DFB_delay_buckets / self.frf

    @property
    def dfb_center_freq(self) -> float:
        """Alias for DFB center frequency [Hz]."""
        return self.frf

    def get_nominal_hc_detuning(
        self,
        I0_mA: float,
        target_xi: float = 0.95,
        use_form_factor: Optional[bool] = None,
    ) -> float:
        """Calculate nominal harmonic cavity detuning in kHz for target xi."""
        I0_A = I0_mA / 1000.0 if I0_mA > 5.0 else I0_mA
        n = self.n_hc
        ratio = self.U0 / self.V_mc
        k_fp_sq = 1.0 / (n**2) - (1.0 / (n**2 - 1.0)) * (ratio**2)
        if k_fp_sq <= 0 or I0_A <= 1e-9 or self.R_hc <= 0 or self.Q_hc <= 0:
            return 35.0  # Default fallback in kHz
        k_fp = np.sqrt(k_fp_sq)
        use_ff = self.use_form_factor if use_form_factor is None else bool(use_form_factor)
        Fn = 0.93 if (use_ff and n >= 3) else 1.0
        if use_ff:
            w_rf = 2.0 * np.pi * self.frf
            tau_max = min(0.45 / self.frf, max(12.0 * self.sigma_t0, 320e-12))
            tau = np.linspace(-tau_max, tau_max, 301)
            dtau = tau[1] - tau[0]
            F1, th1, thn = 0.995, 0.0, 0.0
            for _ in range(18):
                cos_psi = (target_xi * k_fp * self.V_mc) / max(2.0 * Fn * I0_A * self.R_hc, 1e-9)
                psi = 0.0 if cos_psi >= 0.999999 else float(np.arccos(max(0.0, cos_psi)))
                arg_s = np.clip(
                    (self.U0 + 2.0 * (Fn**2) * I0_A * self.R_hc * (np.cos(psi) ** 2)) / max(F1 * self.V_mc, 1.0),
                    -0.999,
                    0.999,
                )
                phis = np.pi - np.arcsin(arg_s) - th1
                Fn_c = Fn * np.exp(-1j * thn)
                E_gain = (
                    self.V_mc * np.sin(phis + w_rf * tau)
                    - 2.0 * I0_A * self.R_hc * np.cos(psi) * np.real(np.exp(1j * psi) * Fn_c * np.exp(1j * n * w_rf * tau))
                    - self.U0
                )
                U_pot = -np.concatenate(([0.0], np.cumsum(0.5 * (E_gain[:-1] + E_gain[1:])) * dtau)) / (self.T0 * self.E0)
                U_pot -= np.min(U_pot)
                rho = np.exp(-np.clip(U_pot / max(self.alpha_c * (self.sigma_delta0**2), 1e-18), 0.0, 100.0))
                rho /= max(float(np.sum(rho) * dtau), 1e-30)
                t_mean = float(np.sum(tau * rho) * dtau)
                tau_c = tau - t_mean
                F1_new = np.sum(rho * np.exp(-1j * w_rf * tau_c)) * dtau
                Fn_new = np.sum(rho * np.exp(-1j * n * w_rf * tau_c)) * dtau
                err = abs(abs(Fn_new) - Fn)
                Fn = 0.5 * float(abs(Fn_new)) + 0.5 * Fn
                F1 = 0.5 * float(abs(F1_new)) + 0.5 * F1
                th1 = float(-np.angle(F1_new))
                thn = float(-np.angle(Fn_new))
                if err < 1e-6:
                    break

        k = (2.0 * Fn * I0_A * self.R_hc) / (target_xi * k_fp * self.V_mc)
        if k <= 1.0:
            return 35.0
        tan_psi = np.sqrt(k**2 - 1.0)
        df_hz = (self.fhc / (2.0 * self.Q_hc)) * tan_psi
        return float(df_hz / 1000.0)

    @property
    def synchrotron_tune_0(self) -> float:
        """Unperturbed linear synchrotron tune Q_s0."""
        # Qs = sqrt(alpha_c * h / (2*pi*E0) * sqrt(V_mc^2 - U0^2))
        cos_phi_s = self.U0 / self.V_mc
        if abs(cos_phi_s) > 1.0:
            return 0.0
        sin_phi_s = np.sqrt(max(0.0, 1.0 - cos_phi_s**2))
        tune_sq = (self.alpha_c * self.h * self.V_mc * sin_phi_s) / (2.0 * np.pi * self.E0)
        return float(np.sqrt(max(0.0, tune_sq)))

    @property
    def synchrotron_freq_0(self) -> float:
        """Unperturbed linear synchrotron frequency f_s0 [Hz]."""
        return self.synchrotron_tune_0 * self.f0

