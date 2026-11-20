function [Vc_mc, Vg_mc_track, Vg_mc_track_0, V_load, V_mc_load_0, PI, DFB] = MC_Cavity_Control(...
    PI, DFB, Vrf_ideal, Vg_mc_track_0, V_mc_load_0, V_load_cpu, TbAng_coef_mc, pattern)
% MC_CAVITY_CONTROL Realistic Main Cavity tracking with independent PI and Direct RF Feedback
%
% This function tracks the main cavity voltage bucket-by-bucket through one turn.
% Supports fully independent operation of:
%   - Proportional-Integral (PI) feedback loop
%   - Direct RF Feedback (DFB) loop
%
% Inputs:
%   PI            : PI controller structure (or [] to disable PI)
%   DFB           : Direct RF Feedback structure (or [] to disable DFB)
%   Vrf_ideal     : Target cavity voltage phasor [V]
%   Vg_mc_track_0 : Generator voltage phasor at start of turn [V]
%   V_mc_load_0   : Beam loading voltage phasor at start of turn [V]
%   V_load_cpu    : Beam-induced voltage for bunches in this turn [V]
%   TbAng_coef_mc : Cavity rotation + decay per bucket
%   pattern       : Fill pattern array (1 = filled bucket, 0 = empty)
%
% Outputs:
%   Vc_mc         : Total cavity voltage phasor array for all buckets in turn [V]
%   Vg_mc_track   : Generator voltage phasor array for all buckets [V]
%   Vg_mc_track_0 : Generator voltage phasor at end of turn [V]
%   V_load        : Beam loading voltage phasor array for all buckets [V]
%   V_mc_load_0   : Beam loading voltage phasor at end of turn [V]
%   PI            : Updated PI structure
%   DFB           : Updated DFB structure

n           = length(pattern);
V_load      = zeros(1, n);
Vg_mc_track = zeros(1, n);
Vc_mc       = zeros(1, n);

use_PI  = ~isempty(PI)  && isfield(PI, 'm')  && (~isfield(PI, 'enable')  || PI.enable);
use_DFB = ~isempty(DFB) && isfield(DFB, 'd') && (~isfield(DFB, 'enable') || DFB.enable);

% Coordinate initial states between PI and DFB to avoid kick mismatch at t=0
if use_PI && use_DFB && PI.piIndex == 0 && isfield(DFB, 'Vg0_main')
    PI.Vg0 = DFB.Vg0_main;
    PI.Ig  = DFB.Ig_main0;
end

j = 1;
for ii = 1:n
    %% 1. PI Feedback / Feedforward Trigger
    if use_PI
        PI.piIndex = PI.piIndex + 1;
        if PI.piIndex > PI.d && mod(PI.piIndex - PI.d, PI.m) == 0
            PI.Vg0 = 2 * pi * PI.RoverQ * PI.Ig_FBid(1);
            PI.Ig_FBid(1) = [];
        end
        Vg0_pi = PI.Vg0;
    elseif ~isempty(PI) && isfield(PI, 'Vg0')
        % PI closed-loop disabled: maintain constant feedforward generator drive
        Vg0_pi = PI.Vg0;
    elseif ~use_DFB && isstruct(DFB) && isfield(DFB, 'Vg0_main')
        % Fallback when PI is empty but DFB has the base generator drive
        Vg0_pi = DFB.Vg0_main;
    else
        Vg0_pi = 0;
    end
    
    %% 2. Direct RF Feedback Trigger
    if use_DFB
        DFB.dfbIndex = DFB.dfbIndex + 1;
        if DFB.dfbIndex > DFB.d && mod(DFB.dfbIndex - DFB.d, DFB.every) == 0
            % Read delayed cavity voltage from FIFO buffer
            vc_delayed = DFB.VcRecord(end);
            
            % Direct negative feedback voltage: vg_drf = -gain * Vc * exp(1j * phase_app)
            % (opposes cavity voltage deviations in STABLE's +Ig -> +Vg phasor convention)
            vg_drf = -DFB.gain * vc_delayed * exp(1i * DFB.phase_shift_app);
            
            % Convert to generator current component (Vg2Ig)
            DFB.Ig_dfb = vg_drf / (2 * DFB.RL * cos(DFB.psi)) * exp(-1i * DFB.psi);
            DFB.Vg0    = 2 * pi * DFB.RoverQ * DFB.Ig_dfb;
            
            if isfield(DFB, 'record_track') && DFB.record_track
                DFB.Ig_track(end+1)     = DFB.Ig_dfb;
                DFB.Vg_drf_track(end+1) = vg_drf;
            end
        end
        Vg0_dfb = DFB.Vg0;
        
        % If PI is not active and PI has no Vg0, include the base feedforward generator drive from DFB
        if ~use_PI && (isempty(PI) || ~isfield(PI, 'Vg0'))
            Vg0_dfb = Vg0_dfb + DFB.Vg0_main;
        end
    else
        Vg0_dfb = 0;
    end
    
    %% 3. Total Generator Kick & Voltage Propagation
    Vg0_total = Vg0_pi + Vg0_dfb;
    Vg_mc_track_0 = Vg_mc_track_0 * TbAng_coef_mc + Vg0_total;
    Vg_mc_track(ii) = Vg_mc_track_0;
    
    %% 4. Beam Loading Voltage Propagation
    if pattern(ii) == 1
        V_load(ii)  = V_mc_load_0 * TbAng_coef_mc;
        V_mc_load_0 = V_load(ii) + V_load_cpu(j);
        j = j + 1;
    else
        V_mc_load_0 = V_mc_load_0 * TbAng_coef_mc;
        V_load(ii)  = V_mc_load_0;
    end
    
    %% 5. Total Cavity Voltage
    Vc_mc(ii) = Vg_mc_track_0 + V_load(ii);
    
    %% 6. Update PI Feedback Sampling & Controller Law
    if use_PI
        piTrackIndex = mod(PI.piIndex, PI.m);
        if piTrackIndex == 0
            piTrackIndex = PI.m;
            PI.Vrf_mc_track(piTrackIndex) = Vc_mc(ii);
            Vrf_mc_mean = mean(PI.Vrf_mc_track);
            
            DI = (Vrf_ideal - Vrf_mc_mean) / PI.RL;
            PI.Integral = PI.Integral + PI.KI * DI;
            PI.DIg      = PI.KP * DI + PI.Integral;
            
            % If DFB is active, PI regulates around the adjusted base current
            if use_DFB && isfield(DFB, 'Ig_main0')
                PI.Ig = DFB.Ig_main0 + PI.DIg;
            else
                PI.Ig = PI.Ig0 + PI.DIg;
            end
            
            PI.Ig_track = [PI.Ig_track, PI.Ig];
            PI.Ig_FBid  = [PI.Ig_FBid,  PI.Ig];
        else
            PI.Vrf_mc_track(piTrackIndex) = Vc_mc(ii);
        end
    end
    
    %% 7. Update Direct RF Feedback Sampling & FIFO Buffer
    if use_DFB
        sample_idx = mod(DFB.dfbIndex - 1, DFB.m) + 1;
        DFB.Vrf_sample_buf(sample_idx) = Vc_mc(ii);
        
        if mod(DFB.dfbIndex, DFB.every) == 0
            valid_len = min(DFB.dfbIndex, DFB.m);
            mean_vc = mean(DFB.Vrf_sample_buf(1:valid_len));
            % Roll FIFO buffer: new value at index 1, oldest dropped from end
            DFB.VcRecord = [mean_vc, DFB.VcRecord(1:end-1)];
        end
    end
end

end
