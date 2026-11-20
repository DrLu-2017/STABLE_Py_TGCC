# Phase Loop Mode 0 Source Index

## Primary implementation location

- `mbtrack2/tracking/rf.py`
  - `class PhaseLoop`

## Key methods and concepts

- `PhaseLoop.__init__(ring, cavity, gain_p=0.0, gain_d=0.0, gain_i=0.0, delay=1, max_correction=None, enabled=True, harmonic=None, mode=None)`
  - Configures the loop gains, delay, control mode, and target cavity.
  - Infers the harmonic number from the cavity if not provided.

- `PhaseLoop.track(bunch)`
  - Called once per turn after the cavity has tracked the bunch.
  - Measures the mean bunch phase, applies delay, computes the control correction, and injects it into the cavity.

- `_measure_phase(bunch_or_beam)`
  - Computes mean `tau` and converts to phase via `m * ring.omega1 * tau_mean`.
  - Handles both single `Bunch` and full `Beam` objects.

- `_compute_correction(phi)`
  - Implements `P`, `PD`, `PI`, and `PID` control laws.
  - Uses backward-Euler derivative and rectangular-rule integral.

- `_apply_correction(dphi)`
  - Directly increments `cavity.theta` for `RFCavity`.
  - For `CavityResonator`, applies correction through `ig_modulation_signal` on feedback loops if available, or falls back to `theta_g`.

- `PhaseLoop.reset()`
  - Resets internal delay buffer, derivative history, integral accumulator, and diagnostics.

- `PhaseLoop.damping_rate_estimate(omega_s)`
  - Returns the approximate damping rate for a PD controller.

## Integration notes

- `PhaseLoop` is a turn-based controller and expects to be called after the RF cavity has updated the bunch.
- The loop can be used with phase-based cavities (`RFCavity`) and beam-loaded cavities (`CavityResonator`).
- The `enabled` flag allows the phase loop to be present in code paths without altering passive behavior.
