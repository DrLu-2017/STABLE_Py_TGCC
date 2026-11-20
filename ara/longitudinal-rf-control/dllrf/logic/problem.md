# DLLRF Problem Statement

The current `mbtrack2` RF feedback infrastructure includes phase-loop and cavity-level feedback models, but it does not yet expose an explicit digital low-level RF (DLLRF) chain in the ARA documentation.

This submodule captures the need to represent the end-to-end digital RF control path used by `CavityResonator`, including:

- ADC sampling at a realistic FPGA clock rate,
- real intermediate-frequency (IF) sampling and digital down-conversion,
- moving-average low-pass filtering,
- FIFO delay in the digital signal path,
- I/Q-domain PI control,
- generator phasor output applied to the cavity.

The implementation is currently contained in `mbtrack2/tracking/rf.py` as the `DLLRFFeedback` class.  This ARA artifact must document the insertion point, the behavior assumptions, and the configuration parameters for the DLLRF model.

## Why this is important

`DLLRFFeedback` is the first feature that makes the RF control model explicitly digital rather than purely phasor-based. Without this documentation, the low-level timing and filtering assumptions are easy to overlook during later refactors or when adding RF phase-noise disturbances.
