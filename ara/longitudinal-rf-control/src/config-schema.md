# Config schema

## Purpose

List the shared and module-specific configuration keys for longitudinal RF control features.

## Shared keys

- `enabled`
- `gain`
- `delay_turns`
- `sign_convention`
- `seed`
- `logging`

## Module-specific examples

### phase_loop

- `observable_type`
- `filter_window`
- `actuator_target`

### dllrf

- `controller_bandwidth`
- `controller_delay`
- `iq_model_enabled`

### rf_phase_noise

- `noise_model`
- `noise_rms`
- `noise_spectrum`
- `injection_target`
