# Validation strategy

## Purpose

Define a layered validation approach shared across all longitudinal RF control modules.

## Validation layers

### Layer 1 — Structural checks

- configuration parsing works;
- feature flags switch behavior cleanly;
- controller/noise states initialize deterministically.

### Layer 2 — Physics smoke tests

- each feature produces a visible and interpretable effect in an isolated case;
- sign and qualitative trend match expectations;
- nominal settings do not create obvious numerical pathologies.

### Layer 3 — Regression checks

- disabling a feature restores baseline behavior;
- existing examples or benchmarks remain within tolerance;
- refactors preserve previously accepted results.

## Shared benchmark philosophy

Each submodule should provide:

- one minimal benchmark;
- one comparison against baseline;
- one artifact in `evidence/` containing plots, metrics, or logs.
