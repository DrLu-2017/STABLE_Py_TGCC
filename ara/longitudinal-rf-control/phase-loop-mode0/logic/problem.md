# Phase Loop Mode 0 Problem Statement

The current `mbtrack2` RF control architecture includes general cavity feedback and a phase loop concept, but it lacks a dedicated ARA artifact documenting the longitudinal mode-0 phase regulator.

This submodule captures the need to make the rigid-bunch (mode-0) phase feedback path explicit, including:

- the measurement of the bunch centroid phase from particle time coordinates (`tau`),
- the conversion from time offset to RF phase via `m * ring.omega1`,
- the control law and delay semantics for PD/PI/PID feedback,
- the actuation path into `RFCavity` (`theta`) and `CavityResonator` (`theta_g` or `ig_modulation_signal`),
- the behavior of the loop under beam loading and with existing cavity feedback objects.

## Why this matters

Without a dedicated mode-0 ARA artifact, the phase control path is easy to conflate with other RF feedback mechanisms:

- `CavityResonator` already supports multiple feedback loops via `cav_res.feedback`, but the mode-0 phase correction is applied differently depending on cavity type.
- The `PhaseLoop` controller is the first explicit mode-0 regulator in the source tree, and its assumptions should be documented before refactor or extension.

## Scope

This artifact should document the `PhaseLoop` class in `mbtrack2/tracking/rf.py`, its integration with both `RFCavity` and `CavityResonator`, and the expected behavior of the loop when enabled or disabled.