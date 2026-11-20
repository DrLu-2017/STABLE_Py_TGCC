# Phase Loop Mode 0 Claims

## Claim 1

The `PhaseLoop` class implements a longitudinal mode-0 phase feedback controller for `mbtrack2`.

## Claim 2

The loop measures the bunch centroid phase from the mean of `tau` and converts it to RF phase using the harmonic number and ring revolution angular frequency.

## Claim 3

The controller supports modes `P`, `PD`, `PI`, and `PID` and applies a user-defined turn delay to maintain causality.

## Claim 4

For `RFCavity`, the phase correction is applied directly to `theta`.
For `CavityResonator`, the correction is injected into the cavity feedback loops via `ig_modulation_signal` if available, otherwise it adjusts `theta_g` directly.

## Claim 5

The loop is configurable and can be disabled cleanly using the `enabled` flag, preserving baseline behavior when off.
