# Source index

## Purpose

This file is the primary navigation map for coding agents and contributors working on longitudinal RF control features in `mbtrack2`.

## Reading order

1. Core longitudinal dynamics update.
2. Beam/bunch observable definitions.
3. RF and cavity state update code.
4. Existing feedback abstractions.
5. Example scripts and configuration patterns.
6. Existing tests or benchmarks.

## Top-level source map

> Replace the placeholders with actual repository paths.

### Core dynamics

- `M:\Project\mbtrack2-dev\mbtrack2\mbtrack2\tracking\particles.py`
- `M:\Project\mbtrack2-dev\mbtrack2\mbtrack2\tracking\element.py`
- `M:\Project\mbtrack2-dev\mbtrack2\mbtrack2\utilities\synchrotron.py`
  
### RF cavity path and feedback class
- `M:\Project\mbtrack2-dev\mbtrack2\mbtrack2\tracking\rf.py`

### Config and examples

- `M:\Project\mbtrack2-dev\mbtrack2\examples\`
- `M:\Project\mbtrack2-dev\mbtrack2\tests\`


## Submodule ownership map

| Submodule | Primary code concern | Likely focus areas |
|---|---|---|
| `phase-loop-mode0` | simple closed-loop phase regulation | observable extraction, delay, actuator insertion |
| `dllrf` | detailed RF controller dynamics | controller state, cavity coupling, timing |
| `rf-phase-noise` | disturbance injection and response | noise source, injection point, reproducible seeds |

## Editing rules for agents

- Modify only one submodule at a time unless a shared refactor is explicitly approved.
- Update trace files before broad refactors.
- Add new config fields before hard-coding parameters.
- If repository inspection reveals an unexpected architecture, update this file first.
