"""
LLRFTuner: LLRF Direct Feedback Crescent-Shaped Stability Window Tuner.

Maps the 2D (K_p vs tau_d) parameter space, extracts the crescent-shaped
stability window, and computes:
1. Chebyshev Center (center of the maximum inscribed sphere / largest margin to boundary)
2. Geometric Centroid (mean center of the stability mask)
Generates publication-quality stability boundary plots with optimal operating points.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Optional, Tuple, Sequence
from pathlib import Path

import numpy as np
from scipy.ndimage import distance_transform_edt

from stable_core.types import MachineParameters, StabilityStatus
from stable_core.classifier import InstabilityClassifier

logger = logging.getLogger(__name__)


@dataclass
class LLRFTuningResult:
    """Results from LLRF feedback parameter space tuning."""
    kp_values: np.ndarray
    delay_values: np.ndarray
    stable_mask: np.ndarray            # 2D boolean mask (len(delay), len(kp))
    nyquist_margins: np.ndarray        # 2D float array
    # Optimal points
    optimal_kp: float                  # Chebyshev center K_p
    optimal_delay: float               # Chebyshev center tau_d
    optimal_margin_score: float        # Distance to boundary / inscribed radius
    centroid_kp: float                 # Geometric centroid K_p
    centroid_delay: float              # Geometric centroid tau_d
    stable_area_fraction: float        # Fraction of parameter space that is stable


class LLRFTuner:
    """
    Tuner for LLRF proportional feedback gain K_p and total loop delay tau_d.
    """

    def __init__(
        self,
        machine: Optional[MachineParameters] = None,
        classifier: Optional[InstabilityClassifier] = None,
    ) -> None:
        self.machine = machine or MachineParameters()
        self.classifier = classifier or InstabilityClassifier(self.machine)

    def scan_feedback_plane(
        self,
        kp_range: Tuple[float, float] = (0.05, 4.0),
        n_kp: int = 50,
        delay_range: Optional[Tuple[float, float]] = None,
        n_delay: int = 50,
        I0_mA: float = 500.0,
        df_khz: float = 40.0,
    ) -> LLRFTuningResult:
        """
        Scan (K_p, tau_d) plane and find optimal operating setpoint.

        Args:
            kp_range: (min_gain, max_gain).
            n_kp: Number of gain steps.
            delay_range: (min_delay, max_delay) in seconds. Defaults to (0.05*T0, 0.95*T0).
            n_delay: Number of delay steps.
            I0_mA: Beam current.
            df_khz: Harmonic cavity detuning.

        Returns:
            LLRFTuningResult with stable mask, Chebyshev center, and centroid.
        """
        if delay_range is None:
            delay_range = (0.05 * self.machine.T0, 2.5 * self.machine.T0)

        kp_arr = np.linspace(kp_range[0], kp_range[1], n_kp)
        delay_arr = np.linspace(delay_range[0], delay_range[1], n_delay)

        stable_mask = np.zeros((n_delay, n_kp), dtype=bool)
        margin_grid = np.zeros((n_delay, n_kp), dtype=float)

        for j, delay in enumerate(delay_arr):
            for i, kp in enumerate(kp_arr):
                # 1. Nyquist loop stability check
                is_nyq, nyq_margin, _ = self.classifier.analyze_llrf_nyquist(kp, delay)

                # 2. Beam stability with feedback
                res = self.classifier.evaluate_semianalytical(
                    I0_mA=I0_mA,
                    df_khz=df_khz,
                    K_p=kp,
                    delay_sec=delay,
                )

                is_stable = is_nyq and res.is_stable and (res.growth_rate <= 0.0)
                stable_mask[j, i] = is_stable
                margin_grid[j, i] = nyq_margin if is_stable else 0.0

        if not np.any(stable_mask):
            logger.warning("No stable points found in scanned LLRF range!")
            return LLRFTuningResult(
                kp_values=kp_arr,
                delay_values=delay_arr,
                stable_mask=stable_mask,
                nyquist_margins=margin_grid,
                optimal_kp=1.0,
                optimal_delay=self.machine.T0 / 2.0,
                optimal_margin_score=0.0,
                centroid_kp=1.0,
                centroid_delay=self.machine.T0 / 2.0,
                stable_area_fraction=0.0,
            )

        # 1. Geometric centroid
        y_indices, x_indices = np.where(stable_mask)
        centroid_j = float(np.mean(y_indices))
        centroid_i = float(np.mean(x_indices))
        centroid_kp = float(np.interp(centroid_i, np.arange(n_kp), kp_arr))
        centroid_delay = float(np.interp(centroid_j, np.arange(n_delay), delay_arr))

        # 2. Chebyshev center using Euclidean Distance Transform (EDT)
        # EDT calculates distance from each stable pixel to the nearest unstable pixel / boundary
        dist_transform = distance_transform_edt(stable_mask)
        max_dist_idx = np.unravel_index(np.argmax(dist_transform), dist_transform.shape)
        chebyshev_j, chebyshev_i = max_dist_idx
        chebyshev_kp = float(kp_arr[chebyshev_i])
        chebyshev_delay = float(delay_arr[chebyshev_j])
        max_inscribed_radius = float(dist_transform[chebyshev_j, chebyshev_i])

        stable_fraction = float(np.sum(stable_mask)) / float(stable_mask.size)

        return LLRFTuningResult(
            kp_values=kp_arr,
            delay_values=delay_arr,
            stable_mask=stable_mask,
            nyquist_margins=margin_grid,
            optimal_kp=chebyshev_kp,
            optimal_delay=chebyshev_delay,
            optimal_margin_score=max_inscribed_radius,
            centroid_kp=centroid_kp,
            centroid_delay=centroid_delay,
            stable_area_fraction=stable_fraction,
        )

    def plot_crescent_window(
        self,
        result: LLRFTuningResult,
        save_path: Optional[str | Path] = None,
        title: str = "LLRF Direct Feedback Stability Window (Crescent Shape)",
    ) -> None:
        """
        Plot 2D crescent stability window, Nyquist margins, and optimal centers.
        """
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

        # Delay in fractions of T0
        delay_frac = result.delay_values / self.machine.T0
        kp = result.kp_values

        # 1. Binary Stability Map
        ax1.pcolormesh(
            kp,
            delay_frac,
            result.stable_mask,
            cmap="RdYlGn",
            shading="auto",
            alpha=0.85,
        )
        ax1.scatter(
            result.optimal_kp,
            result.optimal_delay / self.machine.T0,
            color="blue",
            s=120,
            marker="*",
            edgecolor="black",
            zorder=5,
            label=f"Chebyshev Center (Kp={result.optimal_kp:.2f}, delay={result.optimal_delay/self.machine.T0:.3f} T0)",
        )
        ax1.scatter(
            result.centroid_kp,
            result.centroid_delay / self.machine.T0,
            color="cyan",
            s=80,
            marker="o",
            edgecolor="black",
            zorder=5,
            label=f"Centroid (Kp={result.centroid_kp:.2f}, delay={result.centroid_delay/self.machine.T0:.3f} T0)",
        )
        ax1.set_xlabel("Proportional Gain $K_p$", fontsize=11)
        ax1.set_ylabel(r"Loop Delay $\tau_d / T_0$", fontsize=11)
        ax1.set_title("Crescent Stability Window", fontsize=12, fontweight="bold")
        ax1.grid(True, linestyle="--", alpha=0.5)
        ax1.legend(loc="upper right", fontsize=9)

        # 2. Nyquist Distance Margin Heatmap
        pcm = ax2.pcolormesh(
            kp,
            delay_frac,
            result.nyquist_margins,
            cmap="viridis",
            shading="auto",
        )
        fig.colorbar(pcm, ax=ax2, label="Nyquist Stability Margin")
        ax2.scatter(
            result.optimal_kp,
            result.optimal_delay / self.machine.T0,
            color="red",
            s=120,
            marker="*",
            edgecolor="white",
            zorder=5,
            label="Optimal Working Point",
        )
        ax2.set_xlabel("Proportional Gain $K_p$", fontsize=11)
        ax2.set_ylabel(r"Loop Delay $\tau_d / T_0$", fontsize=11)
        ax2.set_title("Nyquist Margin Landscape", fontsize=12, fontweight="bold")
        ax2.grid(True, linestyle="--", alpha=0.5)
        ax2.legend(loc="upper right", fontsize=9)

        plt.suptitle(title, fontsize=13, fontweight="bold")
        plt.tight_layout()

        if save_path:
            p = Path(save_path)
            p.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(str(p), dpi=200)
            logger.info("Saved LLRF crescent plot to %s", p)
        plt.close(fig)
