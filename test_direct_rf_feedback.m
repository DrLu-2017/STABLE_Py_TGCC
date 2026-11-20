% TEST_DIRECT_RF_FEEDBACK Automated unit and integration test for Direct RF Feedback in STABLE
%
% Tests:
%   1. Standalone DFB_Set initialization and parameter matching.
%   2. Initial steady-state vector matching (Ig_main0 + Ig_dfb0 == Ig_mc_0).
%   3. DFB_calc_properties verification (impedance reduction, gamma, alpha, fs).
%   4. Tracking comparison: PI only vs DFB only vs Combined PI + DFB.
%   5. Backward compatibility of PI_Control with legacy 7-argument call.

clear; clc;
fprintf('=== Starting Direct RF Feedback Test Suite ===\n\n');

%% 1. Setup Machine Parameters (Soleil 2 / HALF test case)
machine_Soleil2;
HALF = machine(C,I0,U0,E0,tau_s,tau_z,sigma_t0,sigma_e0,alpha_c,h,V_mc,n_hc,R_hc,Q_hc,fillrate,fre_shift_hc,Q_mc,R_mc,fre_shift_mc);
HALF.Vb_mc = HALF.qc * HALF.wr_mc * HALF.R_mc / HALF.Q_mc;

rot_decay_coef_mc = 1i * HALF.rot_coef_mc - 1 / (2 * HALF.Q_mc);
TbAng_coef_mc     = exp(rot_decay_coef_mc * HALF.angle_mc);

fprintf('[1/5] Machine initialization complete.\n');
fprintf('      Main Cavity Rs = %.2f kOhm, QL = %.1f, detune angle = %.2f deg\n', ...
    HALF.R_mc/1e3, HALF.Q_mc, HALF.det_angle_mc*180/pi);

%% 2. Test DFB_Set Initialization and Decoupling
DFB_gain        = 0.2;
DFB_phase_shift = 0.0;
DFB_delay       = 800;
DFB_sample_num  = 1;
DFB_every       = 1;

DFB = DFB_Set(HALF, DFB_gain, DFB_phase_shift, DFB_delay, DFB_sample_num, DFB_every);

% Verify DFB structure fields
assert(isfield(DFB, 'gain'), 'DFB must have gain');
assert(isfield(DFB, 'phase_shift'), 'DFB must have phase_shift');
assert(isfield(DFB, 'd'), 'DFB must have delay d');
assert(isfield(DFB, 'm'), 'DFB must have sample_num m');
assert(isfield(DFB, 'every'), 'DFB must have update interval every');
assert(isfield(DFB, 'VcRecord'), 'DFB must have VcRecord buffer');
assert(isfield(DFB, 'Ig_dfb0'), 'DFB must have Ig_dfb0');

% Verify initial steady-state matching: Ig_main0 + Ig_dfb0 must equal HALF.Ig_mc_0
ig_total_0 = DFB.Ig_main0 + DFB.Ig_dfb0;
rel_diff_ig = abs(ig_total_0 - HALF.Ig_mc_0) / abs(HALF.Ig_mc_0);
assert(rel_diff_ig < 1e-12, 'Initial current matching failed: Ig_main0 + Ig_dfb0 != Ig_mc_0');

fprintf('[2/5] DFB_Set successfully initialized.\n');
fprintf('      Ig_mc_0 = %.4f + %.4fi A\n', real(HALF.Ig_mc_0), imag(HALF.Ig_mc_0));
fprintf('      Ig_dfb0 = %.4f + %.4fi A\n', real(DFB.Ig_dfb0), imag(DFB.Ig_dfb0));
fprintf('      Ig_main0= %.4f + %.4fi A\n', real(DFB.Ig_main0), imag(DFB.Ig_main0));
fprintf('      Initial current matching check PASSED (rel error = %.2e)\n', rel_diff_ig);

%% 3. Test DFB_calc_properties
props = DFB_calc_properties(DFB, HALF);

fprintf('[3/5] DFB_calc_properties test:\n');
fprintf('      DFB psi    = %.2f deg\n', props.psi * 180 / pi);
fprintf('      DFB alpha  = %.2f dB\n',  props.alpha);
fprintf('      DFB gamma  = %.4f\n',     props.gamma);
fprintf('      Cavity Rs0 = %.2f kOhm\n', HALF.R_mc / 1e3);
fprintf('      Effective Rs = %.2f kOhm (%.1f%% of uncompensated)\n', ...
    props.Rs / 1e3, props.Rs / HALF.R_mc * 100);
fprintf('      DFB fs     = %.2f kHz\n', props.fs / 1e3);

% Physical sanity: with positive feedback gain, effective Rs should be reduced
assert(props.Rs < HALF.R_mc, 'Effective shunt impedance must be reduced by DFB');
assert(props.gamma > 0, 'DFB gamma must be positive');
assert(props.fs > 0, 'Synchrotron frequency must be positive');
fprintf('      Impedance reduction check PASSED.\n');

%% 4. Multi-turn Cavity Tracking Test: 3 Cases
V_load_cpu_test = double(HALF.Vb_mc * ones(1, sum(pattern==1)));
num_turns_test  = 50;

% Case A: PI only (legacy)
PI_Set;
PI_a = PI;
Vg_mc_track_0_a = HALF.Vg_mc_track_0;
V_mc_load_0_a   = HALF.V_mc_load_0;

for turn = 1:num_turns_test
    [Vc_a, Vg_a, Vg_mc_track_0_a, ~, V_mc_load_0_a, PI_a] = ...
        PI_Control(PI_a, HALF.Vrf_ideal, Vg_mc_track_0_a, V_mc_load_0_a, ...
                   V_load_cpu_test, TbAng_coef_mc, pattern);
end

% Case B: DFB only (PI disabled)
DFB_b = DFB_Set(HALF, 0.2, 0, 800, 1, 1);
Vg_mc_track_0_b = HALF.Vg_mc_track_0;
V_mc_load_0_b   = HALF.V_mc_load_0;

for turn = 1:num_turns_test
    [Vc_b, Vg_b, Vg_mc_track_0_b, ~, V_mc_load_0_b, ~, DFB_b] = ...
        MC_Cavity_Control([], DFB_b, HALF.Vrf_ideal, Vg_mc_track_0_b, V_mc_load_0_b, ...
                          V_load_cpu_test, TbAng_coef_mc, pattern);
end

% Case C: Combined PI + DFB
PI_Set;
PI_c = PI;
DFB_c = DFB_Set(HALF, 0.2, 0, 800, 1, 1);
Vg_mc_track_0_c = HALF.Vg_mc_track_0;
V_mc_load_0_c   = HALF.V_mc_load_0;

for turn = 1:num_turns_test
    [Vc_c, Vg_c, Vg_mc_track_0_c, ~, V_mc_load_0_c, PI_c, DFB_c] = ...
        PI_Control(PI_c, HALF.Vrf_ideal, Vg_mc_track_0_c, V_mc_load_0_c, ...
                   V_load_cpu_test, TbAng_coef_mc, pattern, DFB_c);
end

% Verify voltages remain stable and close to ideal setpoint
err_a = abs(mean(Vc_a) - HALF.Vrf_ideal) / abs(HALF.Vrf_ideal);
err_b = abs(mean(Vc_b) - HALF.Vrf_ideal) / abs(HALF.Vrf_ideal);
err_c = abs(mean(Vc_c) - HALF.Vrf_ideal) / abs(HALF.Vrf_ideal);

fprintf('[4/5] Multi-turn tracking test completed (%d turns):\n', num_turns_test);
fprintf('      Case A (PI only)        rel error to Vrf_ideal: %.4e\n', err_a);
fprintf('      Case B (DFB only)       rel error to Vrf_ideal: %.4e\n', err_b);
fprintf('      Case C (PI + DFB)       rel error to Vrf_ideal: %.4e\n', err_c);

assert(err_a < 0.05, 'PI only voltage error exceeded tolerance');
assert(err_b < 0.05, 'DFB only voltage error exceeded tolerance');
assert(err_c < 0.05, 'Combined PI+DFB voltage error exceeded tolerance');

%% 5. Backward Compatibility Test for PI_Control
% Calling PI_Control with only 6 outputs and 7 inputs (exact legacy signature)
PI_Set;
[Vc_leg, Vg_leg, Vg0_leg, Vl_leg, Vl0_leg, PI_leg] = ...
    PI_Control(PI, HALF.Vrf_ideal, HALF.Vg_mc_track_0, HALF.V_mc_load_0, ...
               V_load_cpu_test, TbAng_coef_mc, pattern);

assert(length(Vc_leg) == length(pattern), 'Legacy call output length mismatch');
fprintf('[5/6] Legacy PI_Control call backward compatibility PASSED.\n');

%% 6. Closed-Loop Perturbation Damping Polarity Test (Negative vs Positive Feedback)
% Measure pure perturbation decay |Vc_pert - Vc_eq| after 5 turns (~0.9 cavity filling times)
DFB_open_eq   = DFB_Set(HALF, 0.0, 0, 704, 1, 1);
DFB_open_pert = DFB_Set(HALF, 0.0, 0, 704, 1, 1);
DFB_neg_eq    = DFB_Set(HALF, 0.3, 0, 704, 1, 1);
DFB_neg_pert  = DFB_Set(HALF, 0.3, 0, 704, 1, 1);

Vg_open_eq   = HALF.Vg_mc_track_0;
Vg_open_pert = HALF.Vg_mc_track_0 + 0.05 * HALF.Vrf_ideal;
Vg_neg_eq    = HALF.Vg_mc_track_0;
Vg_neg_pert  = HALF.Vg_mc_track_0 + 0.05 * HALF.Vrf_ideal;

Vl_open_eq   = HALF.V_mc_load_0;
Vl_open_pert = HALF.V_mc_load_0;
Vl_neg_eq    = HALF.V_mc_load_0;
Vl_neg_pert  = HALF.V_mc_load_0;

for turn = 1:5
    [Vc_open_eq,   ~, Vg_open_eq,   ~, Vl_open_eq,   ~, DFB_open_eq] = ...
        MC_Cavity_Control([], DFB_open_eq,   HALF.Vrf_ideal, Vg_open_eq,   Vl_open_eq,   V_load_cpu_test, TbAng_coef_mc, pattern);
    [Vc_open_pert, ~, Vg_open_pert, ~, Vl_open_pert, ~, DFB_open_pert] = ...
        MC_Cavity_Control([], DFB_open_pert, HALF.Vrf_ideal, Vg_open_pert, Vl_open_pert, V_load_cpu_test, TbAng_coef_mc, pattern);
    [Vc_neg_eq,    ~, Vg_neg_eq,    ~, Vl_neg_eq,    ~, DFB_neg_eq] = ...
        MC_Cavity_Control([], DFB_neg_eq,    HALF.Vrf_ideal, Vg_neg_eq,    Vl_neg_eq,    V_load_cpu_test, TbAng_coef_mc, pattern);
    [Vc_neg_pert,  ~, Vg_neg_pert,  ~, Vl_neg_pert,  ~, DFB_neg_pert] = ...
        MC_Cavity_Control([], DFB_neg_pert,  HALF.Vrf_ideal, Vg_neg_pert,  Vl_neg_pert,  V_load_cpu_test, TbAng_coef_mc, pattern);
end

err_pert_open = abs(mean(Vc_open_pert - Vc_open_eq)) / abs(HALF.Vrf_ideal);
err_pert_neg  = abs(mean(Vc_neg_pert  - Vc_neg_eq))  / abs(HALF.Vrf_ideal);
fprintf('[6/6] Perturbation damping test (5 turns, +5%% initial dVc):\n');
fprintf('      Open-loop (gain=0.0) perturbation remaining: %.4e\n', err_pert_open);
fprintf('      Closed-loop (gain=0.3) perturbation remaining: %.4e\n', err_pert_neg);
assert(err_pert_neg < 0.8 * err_pert_open, ...
    'DFB failed to damp cavity perturbation faster than open-loop (polarity regression!)');
fprintf('      Negative feedback damping polarity check PASSED.\n');

fprintf('\n=== ALL DIRECT RF FEEDBACK TESTS PASSED SUCCESSFULLY! ===\n');
