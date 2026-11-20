# DLLRF Claims

## Claim 1

`DLLRFFeedback` provides an explicit digital low-level RF model for `CavityResonator` within `mbtrack2`.

## Claim 2

The class models multi-rate ADC/DAC timing, real IF sampling, digital down-conversion, low-pass filtering, and FIFO delay before applying generator phasor corrections.

## Claim 3

The DLLRF path is configurable and reusable via the following parameters:

- `f_adc`: ADC sampling frequency,
- `f_if`: intermediate frequency used for DDC,
- `total_delay_s`: total FIFO delay in seconds,
- `lp_taps`: moving-average filter length,
- `kp`, `ki`: PI controller gains.

## Claim 4

When appended to `cav_res.feedback`, `DLLRFFeedback.track()` executes automatically as part of `CavityResonator.track(beam)`, preserving the existing cavity tracking API.

## Claim 5

The documentation in this submodule will make the DLLRF implementation explicit to developers and enable later extension into noise injection or shared controller infrastructure.
