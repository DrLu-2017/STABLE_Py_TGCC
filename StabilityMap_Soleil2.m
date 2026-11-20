% =========================================================================
% StabilityMap_Soleil2.m
%
% 2D stability classification map for Soleil 2
%   X-axis : HHC detuning  (50 to 90 kHz)
%   Y-axis : beam current  (0 to 500 mA)
%   Color  : stable | PTBL | QD-coupling | LCBI | other
%
% Method
% -------
%  For each (Df_hc, I0) point:
%   1. Run short STABLE tracking (Track_num turns).
%   2. Compute centroid oscillation spectrum via FFT.
%   3. Classify instability type from spectral peak frequency and growth
%      rate relative to radiation damping.
%
% Classification criteria (longitudinal, single-cavity approx.)
% ---------------------------------------------------------------
%  STABLE       : max growth rate < damping rate  (no net growth)
%  PTBL         : dominant mode is m=0 (azimuthal, DC current modulation)
%                 -> oscillation at f_s  (synchrotron frequency) or below,
%                    with alternating bunch pattern (anti-phase neighbours)
%  QD coupling  : quadrupole-dipole coupled instability; typically appears
%                 near flat-potential condition; centroid & bunch-length
%                 oscillate simultaneously at f near f_s
%  LCBI         : longitudinal coupled-bunch instability driven by HOM or
%                 HHC; dominant FFT peak at f = k*f_0 � f_s  (k integer)
%                 with growth rate > damping rate
%
% Usage
% -----
%   1. Make sure machine_Soleil2.m, machine.m, PI_Set.m, HOMs_m0.m,
%      Vb_angle_calc.m, TruncatedGaussian.m, VoltageLoadCalc_matlab.m,
%      BinNumCalZ.mexw64, detune_HC_calc.m are all on the MATLAB path.
%   2. Run this script directly.  GPU required.
%   3. After the scan, a comparison figure vs mbtrack2 is generated.
%
% =========================================================================
clc; clear;

%% =========================================================================
%% SCAN GRID  � adjust resolution to your available compute time
%% =========================================================================
Df_scan_kHz = 40:5:120;          % HHC detuning scan [kHz]  (x-axis, widened for n_hc=4)
I0_scan_mA  = 50:50:500;         % beam current scan  [mA]  (y-axis)
% For a quick test, use coarser grids:
% Df_scan_kHz = [40 60 80 100 120];
% I0_scan_mA  = [100 200 300 400 500];

Df_scan = Df_scan_kHz * 1e3;     % [Hz]
I0_scan = I0_scan_mA  * 1e-3;    % [A]

nDf = length(Df_scan);
nI0 = length(I0_scan);

%% =========================================================================
%% FIXED MACHINE PARAMETERS  (Soleil 2)
%% =========================================================================
cspeed   = 299792458;
C        = 354.097;
E0       = 2.75e9;
U0       = 456e3;
h        = 416;
alpha_c  = 1.27e-4;
tau_s    = 7.3e-3;
tau_z    = 7.3e-3;
sigma_e0 = 1.06e-3;
V_mc     = 1.7e6;          % [V]  matching mbtrack2 Vc=1.7e6

n_hc     = 4;              % 4th harmonic (mbtrack2: m_HC=4)
RoQ_hc   = 30;             % R/Q per cavity [Ohm] (mbtrack2: RoQ_HC=30)
Q_hc     = 31e3;           % loaded Q (mbtrack2: Q0_HC=QL_HC=31e3)
Ncav_hc  = 2;              % number of HC cavities (mbtrack2: Ncav_HC=2)
R_hc     = Ncav_hc * Q_hc * RoQ_hc;  % total loaded R = 1.86 MOhm

% MC � matching mbtrack2: Rs=5e6, Q0=35.7e3, QL=6e3, Ncav=4
Rs_mc_per_cav = 5e6;             % shunt impedance per cavity [Ohm]
Q0_mc         = 35.7e3;          % unloaded Q per cavity
Q_mc          = 6e3;             % loaded Q per cavity (mbtrack2: QL=6e3)
Ncav_mc       = 4;               % number of MC cavities (mbtrack2: Ncav_MC=4)
RoQ_mc        = Rs_mc_per_cav / Q0_mc;  % R/Q per cavity ~ 140 Ohm
R_mc          = Ncav_mc * Q_mc * RoQ_mc; % total loaded R ~ 3.36 MOhm
fre_shift_mc  = -35e3;           % [Hz] MC detuning (mbtrack2: MC_det=-35e3)

% Fill pattern � uniform
pattern  = ones(1, h);
fillrate = 1.0;

%% =========================================================================
%% TRACKING PARAMETERS
%% =========================================================================
Par_num    = 1e4;     % macro-particles per bunch  (mbtrack2 uses 1e5; 1e4 for GPU scan)
Track_num  = 10e4;    % total tracking turns (mbtrack2 uses 300e3; 10e4 balances speed/accuracy)
Recor_step = 10;      % record interval [turns]
Recor_num  = Track_num / Recor_step;

% Analysis windows (in recording steps, not turns)
win_transient = round(0.5 * Recor_num);   % discard first 50 % (need >5 damping times to equilibrate)
win_analysis  = win_transient+1 : Recor_num; % use last 50 %

%% =========================================================================
%% INSTABILITY THRESHOLDS
%% =========================================================================
% Growth-rate threshold relative to radiation damping rate (1/tau_z)
GR_threshold_factor = 0.5;   % growth rate > 0.5/tau_z  unstable

% Frequency bands for classification [normalised to f_0 = 1/T0]
% (set after first run; approx values below)
% f_s / f_0 ~ Qs (synchrotron tune) ~ 1e-3 for Soleil 2
Qs_approx   = 1e-3;           % approximate synchrotron tune (update if needed)
f0_approx   = cspeed / C;
fs_approx   = Qs_approx * f0_approx;  % ~ 850 Hz

% Spectral classification windows (fraction of f_rf)
% PTBL  : peak near 0 or f_s, alternating pattern in bunch-by-bunch
% LCBI  : peak at k*f0 � f_s  for any integer k
% QD    : simultaneous growth of mean(Q_std) and mean(Q_mean) near f_s

%% =========================================================================
%% OUTPUT ARRAYS
%% =========================================================================
% Stability code:
%   0 = stable
%   1 = PTBL  (periodic transient bunch lengthening / DC mode)
%   2 = QD coupling
%   3 = LCBI  (coupled-bunch)
%   4 = other / unknown unstable
stability_code = zeros(nI0, nDf);
growth_rate    = zeros(nI0, nDf);   % [1/s]
peak_freq_norm = zeros(nI0, nDf);   % dominant spectral peak / f_s

%% =========================================================================
%% MAIN SCAN LOOP
%% =========================================================================
gd = gpuDevice();
total_points = nI0 * nDf;
pt = 0;

% Pre-generate initial distribution (reused across scan points)
rng(42);
q1 = TruncatedGaussian(1, [-3,3], [Par_num, 1]);
p1 = TruncatedGaussian(1, [-3,3], [Par_num, 1]);

tic;
for iI = 1:nI0
    I0 = I0_scan(iI);

    for iDf = 1:nDf
        fre_shift_hc = Df_scan(iDf);
        pt = pt + 1;

        %% --- Compute natural bunch length for this current ---
        f_rf_est = h * cspeed / C;
        Qs_est   = sqrt(alpha_c / (2*pi) * h * ...
                        sqrt((V_mc/E0)^2 - (U0/E0)^2));
        T0_est   = C / cspeed;
        sigma_t0 = alpha_c * sigma_e0 / (2*pi*Qs_est) * T0_est;

        %% --- Build machine structure ---
        HALF = machine(C, I0, U0, E0, tau_s, tau_z, ...
                       sigma_t0, sigma_e0, alpha_c, ...
                       h, V_mc, ...
                       n_hc, R_hc, Q_hc, ...
                       fillrate, fre_shift_hc, ...
                       Q_mc, R_mc, fre_shift_mc);
        HALF.ShortRange_on = 1;

        %% --- HOMs (none for baseline) ---
        Q_hom_m0 = []; R_hom_m0 = []; fre_hom_m0 = [];
        HOMs_m0;

        %% --- PI feedback ---
        PI_Set;
        PI.KP = 1;
        PI.KI = 1e-4;

        %% --- Bunch charge & beam loading setup ---
        Bun_num = sum(pattern);
        charge  = ones(1,h) .* pattern;
        charge  = charge / sum(charge) * Bun_num;

        HALF.qc     = charge .* pattern * HALF.qc / Par_num;
        HALF.Vb_hc  = HALF.qc * HALF.wr_hc * HALF.R_hc / HALF.Q_hc;
        HALF.Vb_mc  = HALF.qc * HALF.wr_mc * HALF.R_mc / HALF.Q_mc;

        %% --- Initial conditions ---
        Q = gpuArray(single(repmat(q1, 1, Bun_num)));
        P = gpuArray(single(repmat(p1, 1, Bun_num)));
        index_add = gpuArray(single((1:Bun_num) - 1));

        V_hc_load_0 = HALF.V_hc_load_0;
        V_mc_load_0 = HALF.V_mc_load_0;

        %% --- Rotation/decay coefficients ---
        Dq = 0.05;
        rot_decay_coef_hc = 1i*HALF.rot_coef_hc - 1/(2*HALF.Q_hc);
        TbAng_coef_hc     = exp(rot_decay_coef_hc * HALF.angle_hc);
        exp_ang_coef_hc   = -rot_decay_coef_hc * HALF.wr_hc * sigma_t0;

        rot_decay_coef_mc = 1i*HALF.rot_coef_mc - 1/(2*HALF.Q_mc);
        TbAng_coef_mc     = exp(rot_decay_coef_mc * HALF.angle_mc);
        exp_ang_coef_mc   = -rot_decay_coef_mc * HALF.wr_mc * sigma_t0;

        wake_kick_coef = HALF.qc * HALF.kick_coef;

        %% --- Wake function ---
        tau_q = (0:Dq:200)' * sigma_t0;
        Wake_inter_hc = -HALF.wr_hc * R_hc/Q_hc * ...
            exp(-tau_q*HALF.wr_hc/(2*Q_hc)) .* ...
            (cos(tau_q*HALF.wr_hc*HALF.rot_coef_hc) - ...
             HALF.VbImagFactor_hc*sin(tau_q*HALF.wr_hc*HALF.rot_coef_hc));
        Wake_inter_hc(1) = Wake_inter_hc(1)/2;

        Wake_inter_mc = -HALF.wr_mc * R_mc/Q_mc * ...
            exp(-tau_q*HALF.wr_mc/(2*Q_mc)) .* ...
            (cos(tau_q*HALF.wr_mc*HALF.rot_coef_mc) - ...
             HALF.VbImagFactor_mc*sin(tau_q*HALF.wr_mc*HALF.rot_coef_mc));
        Wake_inter_mc(1) = Wake_inter_mc(1)/2;

        Wake_inter = gpuArray(single(Wake_inter_hc + Wake_inter_mc));

        %% --- Recording arrays ---
        record_Q_mean = zeros(Recor_num, Bun_num);
        record_Q_std  = zeros(Recor_num, Bun_num);
        record_P_mean = zeros(Recor_num, Bun_num);
        record_P_std  = zeros(Recor_num, Bun_num);
        Vb_hc_track_record = zeros(1, Track_num);
        record_th = 0;

        HALF.Ig_track      = zeros(1, h/PI.m * Track_num);
        Ig_track_tnum      = 0;
        HALF.Vg_mc_track_0 = HALF.Vg_mc_init;
        Vg_mc              = abs(HALF.Vg_mc_init);
        HALF.rfcoef1_track = HALF.rfcoef1 / HALF.V_mc * Vg_mc;

        %% --- Tracking loop ---
        for i = 1:Track_num
            % drift
            Q = Q + P * HALF.drift_coef;
            Q_min = min(Q);
            Q_new = round((Q - Q_min) * (1/Dq));

            % HHC kick
            exp_angle     = exp(exp_ang_coef_hc * Q);
            exp_angle_sum = gather(sum(exp_angle));
            V_load_cpu    = double(exp_angle_sum .* HALF.Vb_hc(pattern==1));
            [V_load, V_hc_load_0] = VoltageLoadCalc_matlab(V_hc_load_0, V_load_cpu, TbAng_coef_hc, pattern);
            Vb_hc_track_record(i) = mean(V_load);
            V_load_cpu    = V_load(pattern==1) * HALF.kick_coef;
            V_load_gpu    = gpuArray(single(V_load_cpu));
            V_hc_load_kick = V_load_gpu ./ exp_angle;

            % MC kick + PI
            exp_angle     = exp(exp_ang_coef_mc * Q);
            exp_angle_sum = gather(sum(exp_angle));
            V_load_cpu    = double(exp_angle_sum .* HALF.Vb_mc(pattern==1));
            [Vc_mc, Vg_mc_track, HALF.Vg_mc_track_0, V_load, HALF.V_mc_load_0, PI] = ...
                PI_Control(PI, HALF.Vrf_ideal, HALF.Vg_mc_track_0, ...
                           HALF.V_mc_load_0, V_load_cpu, TbAng_coef_mc, pattern);
            V_mc_kick = gpuArray(single(Vc_mc(pattern==1) * HALF.kick_coef)) ./ exp_angle;

            % Short-range wake
            if HALF.ShortRange_on == 1
                binnum = gather(max(max(Q_new))) + 1;
                bin_num_q = sum(BinNumCalZ(binnum, Q_new));
                bin_num_q = reshape(bin_num_q, binnum, Bun_num);
                kick_conv = conv2(bin_num_q, Wake_inter(1:binnum));
                kick_conv = kick_conv(1:binnum,:) .* wake_kick_coef(pattern==1) * min(i/5000, 1);
                Q_new     = Q_new + (1 + index_add * binnum);
                wake_kick = kick_conv(Q_new);
            else
                wake_kick = 0;
            end

            % Radiation damping + quantum excitation
            rad_quan_kick = -HALF.radampcoef * P + ...
                HALF.quanexcoef * gpuArray.randn(Par_num, Bun_num, 'single');

            P = P + rad_quan_kick - HALF.ploss ...
                - real(V_hc_load_kick) + imag(V_hc_load_kick)*HALF.VbImagFactor_hc ...
                - real(V_mc_kick)      + imag(V_mc_kick)     *HALF.VbImagFactor_mc ...
                + wake_kick;

            % Record
            if mod(i, Recor_step) == 0
                record_th = record_th + 1;
                record_Q_mean(record_th,:) = gather(mean(Q));
                record_Q_std(record_th,:)  = gather(std(Q));
                record_P_mean(record_th,:) = gather(mean(P));
                record_P_std(record_th,:)  = gather(std(P));
            end
        end % tracking loop

        %% ----------------------------------------------------------------
        %% POST-PROCESSING: classify instability
        %% ----------------------------------------------------------------
        % --- 1. Growth rate from centroid energy oscillation ---
        %        Use first bunch as representative
        centroid_E = record_P_mean(win_analysis, 1)';
        centroid_Q = record_Q_mean(win_analysis, 1)';
        sigma_BL   = record_Q_std(win_analysis, :) * HALF.sigma_t0;  % [s]

        damping_rate = 1 / HALF.tau_z;   % [1/s]

        % Fit exponential envelope using local maxima (same as GrowthRate_calc.m)
        xx = log(abs(centroid_E) + 1e-30);
        pass_turns = (win_analysis) * Recor_step;
        try
            TF = islocalmax(xx);
            if sum(TF) > 5
                a_fit = polyfit(pass_turns(TF), xx(TF), 1); % fit peaks only
            else
                a_fit = polyfit(pass_turns, xx, 1);         % fallback
            end
            gr_per_turn = a_fit(1);                    % [1/turn]
            gr_si = gr_per_turn / HALF.T0;            % [1/s]
        catch
            gr_si = 0;
        end
        growth_rate(iI, iDf) = gr_si;

        % --- 2. FFT of centroid position (Q_mean) ---
        sig_fft  = centroid_Q - mean(centroid_Q);
        N_fft    = length(sig_fft);
        fft_amp  = abs(fft(sig_fft));
        fft_amp  = fft_amp(1:floor(N_fft/2));
        freq_axis = (0:floor(N_fft/2)-1) / N_fft / (Recor_step * HALF.T0);  % [Hz]

        [~, peak_idx] = max(fft_amp(2:end));  % exclude DC
        peak_idx = peak_idx + 1;
        f_peak_Hz = freq_axis(peak_idx);      % [Hz]

        peak_freq_norm(iI, iDf) = f_peak_Hz / fs_approx;  % normalise to f_s

        % --- 3. Bunch-length oscillation: QD signature ---
        BL_mean_all = mean(sigma_BL, 2);         % average over bunches [s]
        BL_osc_amp  = std(BL_mean_all);           % [s]  oscillation amplitude
        BL_osc_rel  = BL_osc_amp / mean(BL_mean_all);  % relative oscillation

        % --- 4. Bunch-to-bunch pattern: PTBL signature ---
        % Compute bunch-to-bunch centroid variation in last 100 records
        Q_mean_last = mean(record_Q_mean(end-99:end, :), 1);  % [1 x Bun_num]
        std_bb = std(Q_mean_last) * HALF.sigma_t0 * 1e12;     % [ps]

        % ----------------------------------------------------------------
        %  CLASSIFICATION LOGIC
        % ----------------------------------------------------------------
        is_unstable = gr_si > GR_threshold_factor * damping_rate;
        code = 0;   % default: stable

        if is_unstable
            % LCBI: spectral peak well above f_s, likely coupled-bunch mode
            if f_peak_Hz > 3 * fs_approx
                code = 3;   % LCBI

            % QD coupling: both BL and centroid oscillate near f_s,
            %              and BL oscillation is significant
            elseif BL_osc_rel > 0.05 && f_peak_Hz < 3*fs_approx
                code = 2;   % QD coupling

            % PTBL: bunch-by-bunch pattern with low-frequency modulation
            elseif std_bb > 0.5 && f_peak_Hz < 1.5*fs_approx
                code = 1;   % PTBL

            else
                code = 4;   % other / unclassified
            end
        end

        stability_code(iI, iDf) = code;

        %% --- Progress display ---
        fprintf('[%3d/%3d]  I0=%3d mA  Df=%5.1f kHz  -> code=%d  GR=%.2f/tau_z  f_peak=%.0f Hz\n', ...
            pt, total_points, round(I0*1e3), fre_shift_hc*1e-3, code, ...
            gr_si * HALF.tau_z, f_peak_Hz);

    end % Df loop
end % I0 loop

wait(gd);
elapsed = toc;
fprintf('\nScan complete in %.1f min.\n', elapsed/60);

%% =========================================================================
%% SAVE RESULTS
%% =========================================================================
save('StabilityMap_Soleil2_result.mat', ...
     'stability_code', 'growth_rate', 'peak_freq_norm', ...
     'Df_scan_kHz', 'I0_scan_mA', 'HALF');

%% =========================================================================
%% PLOT: Stability Map (STABLE result)
%% =========================================================================
plot_stability_map(Df_scan_kHz, I0_scan_mA, stability_code, ...
    'STABLE � Soleil 2 Stability Map', 1);

%% =========================================================================
%% HELPER FUNCTION � call separately after loading mbtrack2 data
%% =========================================================================
% To overlay mbtrack2 data points, call:
%
%   load('StabilityMap_Soleil2_result.mat');
%   compare_with_mbtrack2(Df_scan_kHz, I0_scan_mA, stability_code, ...
%       mbtrack2_Df, mbtrack2_I0, mbtrack2_code);
%
% where mbtrack2_Df [kHz], mbtrack2_I0 [mA], mbtrack2_code [same numbering]
% are vectors of the data points extracted from the mbtrack2 figure.


% =========================================================================
function plot_stability_map(Df_kHz, I0_mA, code_mat, fig_title, fig_num)
% Plot the 2D stability classification map.
%
%  code: 0=stable, 1=PTBL, 2=QD, 3=LCBI, 4=other

    cmap = [0.20 0.63 0.17;   % 0 stable      � green
            0.99 0.75 0.05;   % 1 PTBL         � amber
            0.85 0.33 0.10;   % 2 QD coupling  � orange-red
            0.69 0.13 0.13;   % 3 LCBI         � dark red
            0.50 0.50 0.50];  % 4 other        � grey

    figure(fig_num); clf;
    set(gcf,'Color','w','Position',[100 100 820 540]);

    imagesc(Df_kHz, I0_mA, code_mat);
    set(gca,'YDir','normal');
    colormap(gca, cmap);
    clim([0 4]);

    cb = colorbar('Ticks', 0.4:1:3.6, ...
                  'TickLabels', {'Stable','PTBL','QD coupling','LCBI','Other'});
    cb.FontSize = 11;

    xlabel('HHC Detuning  \Delta f_{HC}  [kHz]', 'FontSize', 13);
    ylabel('Beam Current  I_0  [mA]',              'FontSize', 13);
    title(fig_title, 'FontSize', 14, 'FontWeight', 'bold');

    set(gca, 'FontSize', 12, 'Box', 'on', 'LineWidth', 1.2, ...
             'XGrid', 'on', 'YGrid', 'on', ...
             'GridColor', [1 1 1]*0.6, 'GridAlpha', 0.4);

    % Optimal detuning curve (analytical, flat-potential condition)
    hold on;
    I0_line  = linspace(min(I0_mA), max(I0_mA), 200) * 1e-3;
    C_ring   = 354.097; h_ring = 416; n_ring = 4;
    U0_ring  = 456e3;   V_mc_r = 1.7e6;
    Q_hc_r   = 31e3;    R_hc_r = 2 * Q_hc_r * 30;  % Ncav=2, RoQ=30
    Df_opt   = arrayfun(@(Ib) detune_HC_calc(Ib,n_ring,C_ring,h_ring, ...
                         U0_ring,V_mc_r,R_hc_r,Q_hc_r), I0_line) * 1e-3;
    plot(Df_opt, I0_line*1e3, 'w--', 'LineWidth', 2.0, ...
         'DisplayName', 'Optimal detuning (flat-pot.)');
    legend('Location','northwest','FontSize',11,'TextColor','w','Color','none');

    drawnow;
end


function compare_with_mbtrack2(Df_kHz, I0_mA, stable_code, ...
                                mb_Df, mb_I0, mb_code)
% Overlay mbtrack2 classified points on top of the STABLE map.
%
% mb_Df   : vector of detuning values [kHz]
% mb_I0   : vector of current values  [mA]
% mb_code : classification code (same scheme as stable_code)

    markers = {'o', 's', '^', 'v', 'd'};   % per code 0..4
    colors  = {'w', 'k', 'k', 'w', 'w'};  % contrast against background

    figure(2); clf;
    set(gcf,'Color','w','Position',[100 100 820 540]);

    % --- Left panel: STABLE ---
    subplot(1,2,1);
    plot_stability_map(Df_kHz, I0_mA, stable_code, 'STABLE', 0);
    % (re-draws into subplot; adjust as needed)

    % --- Right panel: STABLE + mbtrack2 overlay ---
    subplot(1,2,2);
    cmap = [0.20 0.63 0.17; 0.99 0.75 0.05; ...
            0.85 0.33 0.10; 0.69 0.13 0.13; 0.50 0.50 0.50];
    imagesc(Df_kHz, I0_mA, stable_code);
    set(gca,'YDir','normal'); colormap(gca, cmap); clim([0 4]);
    hold on;

    labels = {'Stable','PTBL','QD coupling','LCBI','Other'};
    for k = 0:4
        idx = mb_code == k;
        if any(idx)
            scatter(mb_Df(idx), mb_I0(idx), 80, ...
                    'Marker', markers{k+1}, ...
                    'MarkerEdgeColor', colors{k+1}, ...
                    'MarkerFaceColor', 'none', ...
                    'LineWidth', 2.0, ...
                    'DisplayName', ['mbtrack2: ', labels{k+1}]);
        end
    end

    xlabel('HHC Detuning  \Delta f_{HC}  [kHz]', 'FontSize', 13);
    ylabel('Beam Current  I_0  [mA]',              'FontSize', 13);
    title('STABLE (colour) + mbtrack2 (markers)', 'FontSize', 13, 'FontWeight', 'bold');
    legend('Location','northwest','FontSize',10,'TextColor','w','Color','none');
    set(gca,'FontSize',12,'Box','on','YDir','normal');
    drawnow;
end
