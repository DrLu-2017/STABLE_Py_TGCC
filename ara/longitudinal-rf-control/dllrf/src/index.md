# DLLRF Source Index

## Primary implementation location

- `mbtrack2/tracking/rf.py`
  - `class DLLRFFeedback`

## Key concepts and methods

- `DLLRFFeedback.__init__(ring, cav_res, f_adc=120e6, f_if=10e6, total_delay_s=1e-6, lp_taps=8, kp=0.0, ki=0.0)`
  - Configures the ADC clock, IF frequency, FIFO delay, filter taps, and PI gains.
  - Uses `ring.f1` and `cav_res.m` to compute the RF fundamental frequency.

- `_baseband_to_if(phasor, time_s)`
  - Converts a complex baseband cavity phasor into a real IF sample.

- `_ddc(sample)`
  - Digitally down-converts a real IF sample into I/Q components using a numerically controlled oscillator (NCO).

- `_lowpass(i, q)`
  - Applies a moving-average low-pass filter to the I and Q channels.

- `_compute_iq_correction(i, q)`
  - Computes the complex PI correction in I/Q space.

- `track()`
  - Runs once per turn.
  - Processes all ADC samples within the turn, updates the FIFO delay buffer, and applies the zero-order-hold DAC output to the cavity generator phasor.

- `reset()`
  - Resets internal DLLRF state for replay or repeated experiments.

## Integration notes

- `DLLRFFeedback` is intended to be appended to `cav_res.feedback`.
- `CavityResonator.track(beam)` will call feedback methods on objects in `cav_res.feedback`.
- This implementation is currently embedded in `rf.py`, so the ARA doc is a roadmap for whether it should be refactored into a separate module later.
