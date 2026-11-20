function props = DFB_calc_properties(DFB, HALF, Vc_val, Ig_val, Ig_dfb_val)
% DFB_CALC_PROPERTIES Calculate physical and diagnostic properties of Direct RF Feedback
%
% Implements methods matching mbtrack2's DirectFeedback (direct_rf_feedback.py):
%   - DFB_psi   : Detuning angle with direct RF feedback [rad]
%   - DFB_alpha : Amplitude ratio alpha of DFB [dB]
%   - DFB_gamma : Gain gamma of DFB
%   - DFB_Rs    : Effective cavity shunt impedance [Ohm]
%   - DFB_Vg    : [Vg_main, Vg_drf] generator voltage components [V]
%   - DFB_fs    : Modified synchrotron frequency with DFB [Hz]
%
% References:
%   [1] K. Akai, PRAB 25, 102002 (2022), Eqs. (13), (15).
%   [2] N. Yamamoto et al., Proc. IPAC'23, WEPL161 (2023).

if nargin < 3 || isempty(Vc_val)
    if isfield(HALF, 'Vrf_ideal')
        Vc_val = HALF.Vrf_ideal;
    else
        Vc_val = HALF.V_mc;
    end
end

if nargin < 4 || isempty(Ig_val)
    if isfield(HALF, 'Ig_mc_0')
        Ig_val = HALF.Ig_mc_0;
    else
        Ig_val = 1;
    end
end

if nargin < 5 || isempty(Ig_dfb_val)
    if isfield(DFB, 'Ig_dfb')
        Ig_dfb_val = DFB.Ig_dfb;
    else
        % Compute negative-feedback DFB current from Vc_val
        phase_app = DFB.psi - DFB.phase_shift;
        vg_drf = -DFB.gain * Vc_val * exp(1i * phase_app);
        Ig_dfb_val = vg_drf / (2 * DFB.RL * cos(DFB.psi)) * exp(-1i * DFB.psi);
    end
end

% 1. DFB_psi: Detuning angle with Direct RF Feedback (Fig. 4 of [1])
props.psi = angle(Vc_val) - angle(Ig_val);

% 2. Amplitude ratio factor: fac = |Ig_dfb| / |Ig_main|
% In negative feedback, net Ig = Ig_main + Ig_dfb (where Ig_dfb opposes Vc),
% so the forward reference drive before feedback subtraction is Ig_main = Ig - Ig_dfb.
Ig_main_val = Ig_val - Ig_dfb_val;
fac = abs(Ig_dfb_val) / max(abs(Ig_main_val), eps);

% 3. DFB_alpha: Amplitude ratio alpha in dB (near Eq. (13) of [1])
props.alpha = 20 * log10(max(fac, eps));

% 4. DFB_gamma: Gain gamma of DFB (near Eq. (13) of [1])
if fac < 1
    props.gamma = fac / (1 - fac);
else
    props.gamma = Inf;
end

% 5. DFB_Rs: Effective shunt impedance seen by the beam [Ohm] (Eq. (15) of [1])
props.Rs = HALF.R_mc / (1 + props.gamma * cos(props.psi));

% 6. DFB_Vg: Generator voltage components (Vg_main, Vg_drf)
phase_app = DFB.psi - DFB.phase_shift;
props.Vg_drf = -DFB.gain * Vc_val * exp(1i * phase_app);
if isfield(HALF, 'Vg_mc_init')
    props.Vg_main = HALF.Vg_mc_init - props.Vg_drf;
else
    props.Vg_main = 0;
end

% 7. DFB_fs: Modified synchrotron frequency [Hz]
% Under negative DFB, the lab-frame fixed cavity voltage component seen by a
% dipole beam phase perturbation is Vc_fixed = Vg_main / (1 + gain * exp(1i * phase_app))
% = Vc - V_beam / (1 + gain * exp(1i * phase_app)), suppressing reactive beam loading.
if isfield(HALF, 'fais_mc_whc')
    phi_s = HALF.fais_mc_whc;
else
    phi_s = HALF.fais_nat;
end

if isfield(HALF, 'Vg_mc_init')
    Vc_fixed = props.Vg_main / (1 + DFB.gain * exp(1i * phase_app));
else
    Vc_fixed = Vc_val;
end

th_fixed = angle(Vc_fixed) - angle(Vc_val) + phi_s;
vg_sum   = -abs(Vc_fixed) * cos(th_fixed);

omega_s = 0;
if vg_sum > 0 && isfield(HALF, 'alpha_c') && isfield(HALF, 'w_rf') && isfield(HALF, 'E0') && isfield(HALF, 'T0')
    omega_s = sqrt(HALF.alpha_c * HALF.w_rf * vg_sum / (HALF.E0 * HALF.T0));
end
props.fs = omega_s / (2 * pi);

end
