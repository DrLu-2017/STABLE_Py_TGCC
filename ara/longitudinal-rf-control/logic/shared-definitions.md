# Shared definitions

## Purpose

Provide a common vocabulary so the three submodules do not use the same words with different meanings.

## Suggested terms

- **phase-like observable**: a quantity derived from beam longitudinal state representing phase or arrival-time deviation.
- **actuation point**: the exact place where controller output or disturbance modifies the simulated RF chain.
- **controller state**: internal memory of a feedback or DLLRF block, such as filtered error or delayed output.
- **disturbance source**: an externally defined perturbation entering the beam/RF system, such as RF phase noise.
- **open-loop**: simulation mode with the feedback path disconnected.
- **closed-loop**: simulation mode with controller action applied.
- **baseline**: reference case used for non-regression comparison.

## Naming policy

- Use `phase_loop` for the stage-1 controller.
- Use `dllrf` for low-level RF controller structures.
- Use `rf_phase_noise` for disturbance injection structures.
- Avoid reusing one class name for physically distinct functions.

## Sign convention policy

Every submodule must state clearly:

- the sign of the measured error;
- the sign of the actuator response;
- the coordinate convention used by the underlying longitudinal state.
