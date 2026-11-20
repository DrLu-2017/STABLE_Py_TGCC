function DFB = DFB_Set(HALF, gain, phase_shift, delay, sample_num, every)
% DFB_SET Initialize Direct RF Feedback (DFB) structure for STABLE
%
% Direct RF feedback via generator current (Ig) compensation, independent
% of PI feedback.
%
% Based on mbtrack2 DirectFeedback (direct_rf_feedback.py) and references:
%   [1] K. Akai, "Stability analysis of rf accelerating mode with feedback
%       loops under heavy beam loading in SuperKEKB", PRAB 25, 102002 (2022).
%   [2] N. Yamamoto et al., "Stability survey of a double RF system with RF
%       feedback loops for bunch lengthening in a low-emittance synchrotron
%       ring", Proc. IPAC'23, WEPL161 (2023).
%
% Usage:
%   DFB = DFB_Set(HALF, gain, phase_shift, delay, sample_num, every)
%   or as script: DFB_Set; (uses existing HALF, DFB_gain, etc. in workspace)

if nargin == 0
    % Script mode: read parameters from caller workspace or use defaults
    ws_vars = evalin('caller', 'whos');
    var_names = {ws_vars.name};
    
    if ismember('HALF', var_names)
        HALF = evalin('caller', 'HALF');
    else
        error('DFB_Set: HALF structure must exist in caller workspace.');
    end
    
    if ismember('DFB_gain', var_names)
        gain = evalin('caller', 'DFB_gain');
    else
        gain = 0.2; % Default gain
    end
    
    if ismember('DFB_phase_shift', var_names)
        phase_shift = evalin('caller', 'DFB_phase_shift');
    else
        phase_shift = 0; % Default phase shift in rad
    end
    
    if ismember('DFB_delay', var_names)
        delay = evalin('caller', 'DFB_delay');
    else
        delay = 800; % Default delay in buckets (matches typical STABLE PI delay)
    end
    
    if ismember('DFB_sample_num', var_names)
        sample_num = evalin('caller', 'DFB_sample_num');
    else
        sample_num = 1; % Default sample window in buckets
    end
    
    if ismember('DFB_every', var_names)
        every = evalin('caller', 'DFB_every');
    else
        every = 1; % Default update interval in buckets
    end
else
    if nargin < 2 || isempty(gain), gain = 0.2; end
    if nargin < 3 || isempty(phase_shift), phase_shift = 0; end
    if nargin < 4 || isempty(delay), delay = 800; end
    if nargin < 5 || isempty(sample_num), sample_num = 1; end
    if nargin < 6 || isempty(every), every = 1; end
end

DFB.enable      = true;
DFB.gain        = double(gain);
DFB.phase_shift = double(phase_shift); % User phase shift [rad]
DFB.d           = round(double(delay)); % Latency delay in buckets
DFB.m           = round(double(sample_num)); % Averaging sample window in buckets
DFB.every       = round(double(every)); % Update interval in buckets

if DFB.every < 1
    DFB.every = 1;
end
if DFB.m < 1
    DFB.m = 1;
end
if DFB.d < 1
    DFB.d = 1;
end

% Cavity parameters
DFB.psi            = HALF.det_angle_mc; % Cavity detuning angle [rad]
DFB.phase_shift_app = DFB.psi - DFB.phase_shift; % Applied phase shift: psi - DFB_phase_shift
DFB.RL             = HALF.R_mc;
DFB.QL             = HALF.Q_mc;
DFB.RoverQ         = HALF.R_mc / HALF.Q_mc;

% FIFO delay buffer sizing: ceil(d / every)
DFB.record_size    = ceil(DFB.d / DFB.every);
if DFB.record_size < 1
    DFB.record_size = 1;
end

% Initial steady-state cavity voltage phasor
% (Vrf_ideal or steady cavity phasor)
if isfield(HALF, 'Vrf_ideal') && ~isempty(HALF.Vrf_ideal)
    vc_init = HALF.Vrf_ideal;
elseif isfield(HALF, 'V_mc') && ~isempty(HALF.V_mc)
    vc_init = -HALF.V_mc * exp(-1i*(pi/2 - HALF.fais_mc_whc));
else
    vc_init = 0;
end

% Initial direct RF feedback voltage and current components
% Negative feedback polarity: in STABLE's cavity phasor convention (machine.m),
% Vg = +2 * RL * cos(psi) * exp(+1i * psi) * Ig and Vc = Vg + V_load, so opposing
% cavity voltage perturbations dVc requires a minus sign (-gain).
vg_drf_init        = -DFB.gain * vc_init * exp(1i * DFB.phase_shift_app);

% Vg2Ig conversion:
% In STABLE / accelerator cavity model:
% Ig = Vg / (2 * RL * cos(psi)) * exp(-1i * psi)
% which is Vg * (1 - 1i*tan(psi)) / (2 * RL)
DFB.Ig_dfb0        = vg_drf_init / (2 * DFB.RL * cos(DFB.psi)) * exp(-1i * DFB.psi);
DFB.Vg0_dfb        = 2 * pi * DFB.RoverQ * DFB.Ig_dfb0;

% Pre-adjust base generator current to match cavity at t=0 (avoid transient shock)
% Ig_main = Ig0 - Ig_dfb0
if isfield(HALF, 'Ig_mc_0')
    DFB.Ig_main0   = HALF.Ig_mc_0 - DFB.Ig_dfb0;
    DFB.Vg0_main   = 2 * pi * DFB.RoverQ * DFB.Ig_main0;
else
    DFB.Ig_main0   = 0;
    DFB.Vg0_main   = 0;
end

% Buffers for delay and averaging
DFB.VcRecord       = ones(1, DFB.record_size) * vc_init;
DFB.vc_previous    = ones(1, DFB.m) * vc_init;
DFB.Vrf_sample_buf = zeros(1, DFB.m);

% Runtime tracking states
DFB.dfbIndex       = 0;
DFB.Ig_dfb         = DFB.Ig_dfb0;
DFB.Vg0            = DFB.Vg0_dfb;
DFB.Ig_track       = [];
DFB.Vg_drf_track   = [];
DFB.record_track   = false; % Set to true only if bucket-by-bucket DFB history is required

% In script mode, assign DFB into caller workspace
if nargout == 0
    assignin('caller', 'DFB', DFB);
end
end
