# Integration Roadmap

Planned integration steps and milestones for the RF control system.# Integration roadmap

## Strategy

The implementation strategy is incremental. Features are added in an order that reduces ambiguity and supports progressive validation.

## Stage 1 — Phase loop mode 0

### Goal

Implement the simplest closed-loop longitudinal phase controller acting on the mode 0 component.

### Why first

This stage establishes:

- the measurement path for a phase-like longitudinal observable;
- the basic controller scheduling pattern;
- the RF-related actuation interface;
- the initial validation workflow for open-loop versus closed-loop comparison.

### Exit criteria

- feature is configurable;
- feature can be enabled or disabled cleanly;
- one benchmark demonstrates intended damping or regulation behavior;
- baseline cases remain unchanged when disabled.

## Stage 2 — DLLRF

### Goal

Introduce a more realistic low-level RF control structure on top of, or alongside, the stage-1 control path.

### Questions to resolve

- which controller states belong in DLLRF rather than in a simplified phase loop;
- whether cavity dynamics and control dynamics should be separated explicitly;
- whether delay and filtering already introduced in stage 1 should be refactored into shared infrastructure.

### Exit criteria

- DLLRF signal path is explicit;
- controller states and delays are documented;
- one benchmark demonstrates meaningful closed-loop behavior;
- interactions with phase loop are documented.

## Stage 3 — RF phase noise

### Goal

Introduce RF phase noise as a disturbance source and evaluate its propagation through the RF/beam system.

### Questions to resolve

- where noise is injected;
- whether noise is applied to reference phase, cavity phase, or drive phase;
- how to compare noisy and noiseless cases reproducibly.

### Exit criteria

- noise model is parameterized;
- injection point is explicit;
- at least one sensitivity or response benchmark is available;
- interactions with phase loop and DLLRF are documented.

## Integration policy

A later stage may refactor earlier code only if:

- behavior is preserved under a documented regression check;
- the reason for refactor is recorded in `trace/master-decisions.yaml`;
- submodule-specific trace files record local consequences.
