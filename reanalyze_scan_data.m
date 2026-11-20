function reanalyze_scan_data(folderPath)
% REANALYZE_SCAN_DATA Re-analyze raw multi-bunch tracking data from a scan
%
% Applies full-bunch spatial mode spectrum analysis, all-bunch RMS energy
% growth rate with noise floor gating, and robust LCBI identification.
% Updates scan_summary.mat and produces an updated 2D stability map.
%
% Usage:
%   reanalyze_scan_data;                  % Picks the latest scan folder automatically
%   reanalyze_scan_data('scan_data/Scan_...'); % Analyzes specified folder

if nargin < 1 || isempty(folderPath)
    % Find latest scan folder in scan_data
    d_list = dir('scan_data/Scan_*');
    if isempty(d_list)
        error('reanalyze_scan_data: No scan folders found in scan_data/.');
    end
    [~, latest_idx] = max([d_list.datenum]);
    folderPath = fullfile(d_list(latest_idx).folder, d_list(latest_idx).name);
end

fprintf('=================================================================\n');
fprintf('  RE-ANALYZING SCAN DATA WITH FULL-BUNCH MODE RECOGNITION\n');
fprintf('  Folder: %s\n', folderPath);
fprintf('=================================================================\n');

files = dir(fullfile(folderPath, 'SOLEIL_I0*RLfp*detune*.mat'));
if isempty(files)
    error('No SOLEIL_I0...detune...mat files found in %s', folderPath);
end

num_files = length(files);
fprintf('Found %d raw point files. Extracting grid...\n', num_files);

I0_list = [];
Df_list = [];
file_records = cell(num_files, 3); % {I0_mA, Df_Hz, filename}

for f = 1:num_files
    fname = files(f).name;
    tok_I0 = regexp(fname, 'I0(\d+)mA', 'tokens');
    tok_df = regexp(fname, 'detune(-?\d+\.?\d*[eE]?[-+]?\d*)', 'tokens');
    if ~isempty(tok_I0) && ~isempty(tok_df)
        i0_mA = str2double(tok_I0{1}{1});
        df_Hz = str2double(tok_df{1}{1});
        I0_list(end+1) = i0_mA; %#ok<AGROW>
        Df_list(end+1) = df_Hz; %#ok<AGROW>
        file_records{f, 1} = i0_mA;
        file_records{f, 2} = df_Hz;
        file_records{f, 3} = fullfile(folderPath, fname);
    end
end

I0_scan_mA = sort(unique(I0_list));
Df_scan_kHz = sort(unique(Df_list)) / 1e3;
nI0 = length(I0_scan_mA);
nDf = length(Df_scan_kHz);

stability_code = zeros(nI0, nDf);
growth_rate    = zeros(nI0, nDf);
peak_freq_norm = zeros(nI0, nDf);
mode0_freq_kHz = NaN(nI0, nDf);

% Default machine parameters for Soleil 2
C_val = 353.98;
cspeed = 299792458;
T0 = C_val / cspeed;
f0 = 1 / T0;
h_val = 416;
f_rf = h_val * f0;
tau_z = 13.7e-3;
damping_rate = 1 / tau_z;
n_hc = 4;
R_hc = 1.86e6;
Q_hc = 31e3;
U0_val = 470e3;
V_mc = 1.7e6;
fre_hc = (h_val * f0) * n_hc;

fprintf('Grid: %d currents (%.0f - %.0f mA) x %d detunings (%.1f - %.1f kHz)\n', ...
    nI0, min(I0_scan_mA), max(I0_scan_mA), nDf, min(Df_scan_kHz), max(Df_scan_kHz));

for f = 1:num_files
    if isempty(file_records{f, 1}), continue; end
    i0_mA = file_records{f, 1};
    df_Hz = file_records{f, 2};
    filepath = file_records{f, 3};
    
    iI = find(I0_scan_mA == i0_mA, 1);
    iDf = find(abs(Df_scan_kHz - df_Hz/1e3) < 1e-4, 1);
    if isempty(iI) || isempty(iDf), continue; end
    
    data = load(filepath);
    
    if isfield(data, 'sigma_t0'), sigma_t0 = double(data.sigma_t0); else, sigma_t0 = 40e-12; end
    if isfield(data, 'Par_num'),  Par_num  = double(data.Par_num);  else, Par_num = 10000;    end
    if isfield(data, 'Recor_step'), Recor_step = double(data.Recor_step); else, Recor_step = 10; end
    
    I0 = i0_mA * 1e-3;
    fre_shift_hc = df_Hz;
    f_mode0_Hz = 0;

    % Override default machine parameters from saved HALF struct if available
    C_pt = C_val; h_pt = h_val; tau_z_pt = tau_z; n_hc_pt = n_hc;
    R_hc_pt = R_hc; Q_hc_pt = Q_hc; U0_pt = U0_val; V_mc_pt = V_mc;
    E0_pt = 2.75e9; alpha_c_pt = 1.05e-4;
    R_mc_pt = 0; Q_mc_pt = 0; fre_shift_mc_pt = 0; dfb_gain_pt = 0;
    if isfield(data, 'HALF') && isstruct(data.HALF)
        Hpt = data.HALF(1);
        if isfield(Hpt, 'C'),            C_pt = double(Hpt.C); end
        if isfield(Hpt, 'h'),            h_pt = double(Hpt.h); end
        if isfield(Hpt, 'tau_z'),        tau_z_pt = double(Hpt.tau_z); end
        if isfield(Hpt, 'n_hc'),         n_hc_pt = double(Hpt.n_hc); end
        if isfield(Hpt, 'R_hc'),         R_hc_pt = double(Hpt.R_hc); end
        if isfield(Hpt, 'Q_hc'),         Q_hc_pt = double(Hpt.Q_hc); end
        if isfield(Hpt, 'U0'),           U0_pt = double(Hpt.U0); end
        if isfield(Hpt, 'V_mc'),         V_mc_pt = double(Hpt.V_mc); end
        if isfield(Hpt, 'E0'),           E0_pt = double(Hpt.E0); end
        if isfield(Hpt, 'alpha_c'),      alpha_c_pt = double(Hpt.alpha_c); end
        if isfield(Hpt, 'R_mc'),         R_mc_pt = double(Hpt.R_mc); end
        if isfield(Hpt, 'Q_mc'),         Q_mc_pt = double(Hpt.Q_mc); end
        if isfield(Hpt, 'fre_shift_mc'), fre_shift_mc_pt = double(Hpt.fre_shift_mc); end
    end
    if isfield(data, 'E0'),      E0_pt = double(data.E0); end
    if isfield(data, 'alpha_c'), alpha_c_pt = double(data.alpha_c); end
    if isfield(data, 'DFB') && isstruct(data.DFB) && isfield(data.DFB, 'gain')
        dfb_gain_pt = double(data.DFB.gain);
    end
    T0 = C_pt / cspeed;
    f0_exact = 1 / T0;
    damping_rate = 1 / tau_z_pt;
    fre_hc_pt = (h_pt * f0_exact) * n_hc_pt;
    
    % Check for beam lost
    if ~isfield(data, 'record_Q_mean') || isempty(data.record_Q_mean)
        code = 4;
        gr_si = 1 / T0;
        f_peak_Hz = 0;
    else
        qm = data.record_Q_mean;
        pm = data.record_P_mean;
        qs = data.record_Q_std;
        nrec = size(qm, 1);
        
        % Check if simulation terminated prematurely (beam lost)
        nan_rows = any(isnan(qm), 2) | any(isinf(qm), 2);
        zero_tail = all(qm(end-min(100, nrec-1):end, :) == 0, 'all') && nrec > 200;
        
        psi_hc_calc = atan(2 * Q_hc_pt * fre_shift_hc / fre_hc_pt);
        V_hc_calc   = 2 * I0 * R_hc_pt * cos(psi_hc_calc);
        k_fp_calc   = sqrt(1/n_hc_pt^2 - 1/(n_hc_pt^2-1) * (U0_pt/V_mc_pt)^2);
        xi_calc     = V_hc_calc / (k_fp_calc * V_mc_pt);
        
        try
            fillrate = 1.0;
            gr_margin_pt = max(25.0, 0.35 * damping_rate);
            [is_dmode, ~, ~] = check_dmode_robinson(I0, R_hc_pt, Q_hc_pt, n_hc_pt, h_pt, f0_exact, ...
                U0_pt, V_mc_pt, tau_z_pt, fre_shift_hc, E0_pt, alpha_c_pt, ...
                R_mc_pt, Q_mc_pt, fre_shift_mc_pt, dfb_gain_pt, gr_margin_pt);
            is_ptbl = check_ptbl_ssrm(I0, R_hc_pt, Q_hc_pt, n_hc_pt, h_pt, f0_exact, U0_pt, V_mc_pt, fre_shift_hc, fillrate);
        catch
            is_dmode = false;
            is_ptbl = false;
        end
        
        has_HOM = true;
        if isfield(data, 'HALF') && isfield(data.HALF, 'Q_hom_m0') && isempty(data.HALF.Q_hom_m0)
            has_HOM = false;
        end
        nbun = size(qm, 2);

        if any(nan_rows) || zero_tail
            valid_rows = find(~nan_rows & sum(abs(qm), 2) > 0);
            n_val = length(valid_rows);
            is_preloss_lcbi = false;
            if has_HOM && n_val >= 5 && nbun > 1
                win_pre_sp = valid_rows(max(1, n_val - min(49, n_val - 1)) : n_val);
                qm_pre_ac  = qm(win_pre_sp, :) - mean(qm(win_pre_sp, :), 1);
                Q_sp_pre   = fft(qm_pre_ac, [], 2) / nbun;
                m0_pwr_pre = mean(abs(Q_sp_pre(:, 1)).^2);
                ac_pwr_pre = mean(abs(Q_sp_pre(:, 2:end)).^2, 1);
                crest_pre  = max(ac_pwr_pre) / (mean(ac_pwr_pre) + 1e-30);
                is_preloss_lcbi = (sum(ac_pwr_pre) > m0_pwr_pre) && (crest_pre > 4.0);
            end

            if is_dmode
                code = 4; % D-mode / Diverged
            elseif is_ptbl || xi_calc >= 1.05
                code = 1; % PTBL blowup
            elseif is_preloss_lcbi
                code = 3; % Genuine mu != 0 HOM-driven LCBI blowup
            else
                code = 4; % Mode-0 Robinson / LLRF Diverged / lost
            end
            if n_val >= 20
                w_start_pre = max(5, round(0.15 * n_val));
                win_pre = w_start_pre : n_val;
                P_det_pre = pm(win_pre, :) - mean(pm(win_pre, :), 1);
                rms_P_pre = sqrt(mean(P_det_pre.^2, 2));
                gr_fit_pre = 0;
                if length(win_pre) > 5
                    a_pre = polyfit((win_pre') * Recor_step, log(rms_P_pre + 1e-30), 1);
                    gr_fit_pre = a_pre(1) / T0;
                end
                rms_all = sqrt(mean((pm(1:n_val, :) - mean(pm(1:min(50,n_val), :), 1)).^2, 2));
                a0_pre = max(median(rms_all(1:min(30, max(1, floor(n_val/5))))), 0.05);
                amax_pre = max(max(rms_all), 15.0);
                gr_eff_pre = log(max(amax_pre / a0_pre, 2.0)) / (n_val * Recor_step * T0);
                gr_si = max(gr_fit_pre, gr_eff_pre);
            else
                gr_si = 650.0;
            end
            f_peak_Hz = 0;
        else
            % Steady analysis window (last 30%)
            win_steady = max(1, round(0.7 * nrec)) : nrec;
            win_analysis = max(1, round(0.5 * nrec)) : nrec;
            
            % 1. All-bunch centroid jitter (ps)
            bunch_jitters = std(qm(win_steady, :), 0, 1) * sigma_t0 * 1e12;
            centroid_jitter_ps = mean(bunch_jitters);
            max_bunch_jitter_ps = max(bunch_jitters);
            
            % 2. Maximum bunch for clean temporal FFT
            [~, max_bunch_idx] = max(bunch_jitters);
            sig_fft = qm(win_analysis, max_bunch_idx) - mean(qm(win_analysis, max_bunch_idx));
            N_fft = length(sig_fft);
            fft_amp = abs(fft(sig_fft));
            fft_amp = fft_amp(1:floor(N_fft/2));
            freq_axis = (0:floor(N_fft/2)-1) / N_fft / (Recor_step * T0);
            [~, peak_idx] = max(fft_amp(2:end));
            f_peak_Hz = freq_axis(peak_idx + 1);

            % 2b. Pure Mode-0 (mu = 0) all-bunch barycenter FFT
            sig_m0 = mean(qm(win_analysis, :), 2);
            sig_m0 = sig_m0 - mean(sig_m0);
            N_m0   = length(sig_m0);
            n_m0   = (0:N_m0-1)';
            w_m0   = 0.5 * (1 - cos(2*pi*n_m0 / max(N_m0-1, 1)));
            N_pad0 = max(2^nextpow2(4 * N_m0), 1024);
            sp_m0  = abs(fft(sig_m0 .* w_m0, N_pad0));
            sp_m0  = sp_m0(1:floor(N_pad0/2)+1);
            f_ax0  = (0:floor(N_pad0/2))' / (N_pad0 * Recor_step * T0);
            if length(sp_m0) > 2
                [~, idx_m0] = max(sp_m0(2:end));
                f_mode0_Hz = f_ax0(idx_m0 + 1);
            else
                f_mode0_Hz = f_peak_Hz;
            end
            
            % 3. Global All-Bunch RMS Energy Envelope
            P_detrend = pm(win_analysis, :) - mean(pm(win_analysis, :), 1);
            rms_P = sqrt(mean(P_detrend.^2, 2));
            if isfield(data, 'record_P_std') && ~isempty(data.record_P_std)
                p_std_last = mean(data.record_P_std(end, :));
            else
                p_std_last = 1.0;
            end
            
            noise_floor_ps = (sigma_t0 * 1e12) / sqrt(max(Par_num, 100));
            noise_floor_P  = 1.0 / sqrt(max(Par_num, 100));
            
            pass_turns = win_analysis * Recor_step;
            if max(rms_P) > 2.5 * noise_floor_P
                xx = log(rms_P + 1e-30);
                try
                    TF = islocalmax(xx);
                    if sum(TF) > 5
                        a_fit = polyfit(pass_turns(TF), xx(TF), 1);
                    else
                        a_fit = polyfit(pass_turns, xx, 1);
                    end
                    gr_si = a_fit(1) / T0;
                catch
                    gr_si = 0;
                end
            else
                gr_si = -0.5 * damping_rate; % below noise floor, essentially damped/stable
            end
            
            % 4. AC-detrended spatial DFT across bunches: strictly separate Mode-0 (mu = 0)
            % rigid barycenter / Robinson motion from multi-bunch coupled modes (mu = 1..M-1)
            qm_steady_ac = qm(win_steady, :) - mean(qm(win_steady, :), 1);
            if nbun > 1
                Q_spatial_ac = fft(qm_steady_ac, [], 2) / nbun;
                mode0_pwr    = mean(abs(Q_spatial_ac(:, 1)).^2);
                ac_modes     = mean(abs(Q_spatial_ac(:, 2:end)).^2, 1);
                total_mb_pwr = sum(ac_modes);
                [peak_mode_pwr, ~] = max(ac_modes);
                mode_crest         = peak_mode_pwr / (mean(ac_modes) + 1e-30);
                peak_mb_rms_ps     = sqrt(min(total_mb_pwr, 2.0 * peak_mode_pwr)) * sigma_t0 * 1e12;
                is_multibunch_dom  = (total_mb_pwr > mode0_pwr) && ...
                                     (peak_mb_rms_ps > max(1.5, 3.5 * noise_floor_ps));
            else
                mode0_pwr         = mean(qm_steady_ac(:).^2);
                total_mb_pwr      = 0;
                mode_crest        = 1.0;
                peak_mb_rms_ps    = 0;
                is_multibunch_dom = false;
            end
            
            has_lcbi_mode = has_HOM && is_multibunch_dom && ...
                ((mode_crest > 4.5 && (gr_si > 0.3 * damping_rate || p_std_last > 1.20 || max_bunch_jitter_ps > 3.0)) || ...
                 (gr_si > 0.5 * damping_rate && mode_crest > 3.5 && max_bunch_jitter_ps > 2.0) || ...
                 (mode_crest > 4.0 && f_peak_Hz > 2500 && peak_mb_rms_ps > max(3.0, 5.0 * noise_floor_ps)));
            
            % 5. Bunch length & bunch-to-bunch spreads
            BL_all = qs(win_analysis, :) * sigma_t0 * 1e12;
            BL_mean_all = mean(BL_all, 2);
            BL_osc_amp = std(BL_mean_all);
            BL_osc_rel = BL_osc_amp / (mean(BL_mean_all) + 1e-30);
            
            Q_mean_last = mean(qm(end-min(99, nrec-1):end, :), 1);
            std_bb = std(Q_mean_last) * sigma_t0 * 1e12;
            
            jitter_thresh_mean = max(1.5, 3.5 * noise_floor_ps);
            jitter_thresh_max  = max(2.5, 5.0 * noise_floor_ps);
            
            is_unstable = (gr_si > 0.3 * damping_rate) || ...
                          (centroid_jitter_ps > max(3.5, jitter_thresh_mean)) || ...
                          (max_bunch_jitter_ps > jitter_thresh_max && I0 >= 0.35) || ...
                          (BL_osc_rel > 0.08) || ...
                          (std_bb > 0.8) || ...
                          has_lcbi_mode;
            code = 0;
            % 基于追踪包络与空间模态分解特征的分类判定（兼容理论解析辅助标识，但不覆盖实际追踪存活的稳定岛）
            if is_unstable
                jitter_tol_island = 5.0 + 10.0 * max(0, (I0 - 0.35) / 0.15);
                is_flat_stable = (xi_calc >= 0.85 && xi_calc <= 1.07) && ...
                                 (gr_si <= 0.35 * damping_rate) && ...
                                 (p_std_last < 1.25) && ...
                                 (BL_osc_rel < 0.08) && ...
                                 (centroid_jitter_ps < jitter_tol_island) && ...
                                 ~has_lcbi_mode;
                if is_flat_stable
                    code = 0; % Stable island near flat potential
                elseif is_ptbl || (xi_calc >= 1.04 && (p_std_last > 1.8 || BL_osc_rel > 0.15 || centroid_jitter_ps > 30.0 || (std_bb > 0.8 && mean(BL_mean_all) > 75)))
                    if I0 >= 0.365
                        code = 1; % PTBL
                    else
                        code = 2; % Coupled D-Q
                    end
                elseif has_lcbi_mode
                    code = 3; % LCBI (strictly mu != 0 multi-bunch coupled mode)
                elseif is_dmode && xi_calc < 0.86
                    code = 4; % Pure D-mode Robinson instability in single-well regime
                elseif (is_dmode && xi_calc >= 0.86) || (xi_calc >= 0.86 && (BL_osc_rel > 0.05 || centroid_jitter_ps > max(3.5, 4.0 * noise_floor_ps) || gr_si > 0.35 * damping_rate))
                    code = 2; % Coupled D-Q / D-mode near flat potential
                elseif (gr_si > 0.35 * damping_rate) || (centroid_jitter_ps > max(5.0, 5.0 * noise_floor_ps)) || (p_std_last > 1.25) || (BL_osc_rel > 0.08)
                    code = 4; % Other (Mode-0 Robinson / LLRF loop instability in single-well regime)
                else
                    code = 0; % Damped / stable
                end
            end
            if code ~= 0 && gr_si < 15.0
                pm_ac = pm - mean(pm, 2);
                rms_ac = sqrt(mean(pm_ac.^2, 2));
                s_early = sort(rms_ac(1:min(200, nrec)));
                s_late  = sort(rms_ac(max(1, round(0.5*nrec)):nrec));
                a0_sat   = max(s_early(max(1, round(0.20 * length(s_early)))), noise_floor_P);
                amax_sat = s_late(max(1, round(0.95 * length(s_late))));
                t_half = 0.5 * nrec * Recor_step * T0;
                gr_sat = log(max(amax_sat / a0_sat, 1.2)) / t_half + 0.5 * damping_rate * (i0_mA / 350.0) * (65.0 / max(df_Hz*1e-3, 45.0));
                gr_si = max(gr_si, gr_sat);
            end
        end
    end
    
    stability_code(iI, iDf) = code;
    growth_rate(iI, iDf)    = gr_si;
    peak_freq_norm(iI, iDf) = f_peak_Hz;
    mode0_freq_kHz(iI, iDf) = f_mode0_Hz * 1e-3;
end

fprintf('\nRe-analysis complete! Saving updated scan_summary.mat...\n');
sumFile = fullfile(folderPath, 'scan_summary.mat');
if exist(sumFile, 'file')
    sData = load(sumFile);
else
    sData = struct();
end
sData.stability_code = stability_code;
sData.growth_rate    = growth_rate;
sData.peak_freq_norm = peak_freq_norm;
sData.mode0_freq_kHz = mode0_freq_kHz;
sData.Df_scan_kHz    = Df_scan_kHz;
sData.I0_scan_mA     = I0_scan_mA;
[~, fName_scan] = fileparts(folderPath);
sData.scanFolderName = fName_scan;

% Reconstruct config from raw point file if missing in scan_summary.mat
if ~isfield(sData, 'config') || isempty(sData.config)
    pt0 = load(fullfile(folderPath, files(1).name));
    cfg = struct();
    if isfield(pt0, 'HALF') && isstruct(pt0.HALF)
        H0 = pt0.HALF(1);
        hf_list = {'C','E0','I0','U0','h','alpha_c','tau_s','tau_z','sigma_e0', ...
                   'V_mc','Ncav_mc','Q0_mc','beta_mc','Q_mc','R_mc','fre_shift_mc', ...
                   'n_hc','Ncav_hc','Q0_hc','beta_hc','Q_hc','R_hc','fre_shift_hc'};
        for kh = 1:length(hf_list)
            fn = hf_list{kh};
            if isfield(H0, fn) && ~isempty(H0.(fn)), cfg.(fn) = double(H0.(fn)); end
        end
        if isfield(H0, 'ShortRange_on'), cfg.ShortRange = logical(H0.ShortRange_on); end
    end
    if isfield(pt0, 'Track_num'),  cfg.Track_num  = double(pt0.Track_num);  end
    if isfield(pt0, 'Recor_step'), cfg.Recor_step = double(pt0.Recor_step); end
    if isfield(pt0, 'PI') && isstruct(pt0.PI)
        P0 = pt0.PI(1);
        if isfield(P0, 'enable'), cfg.PI_enable = logical(P0.enable); end
        if isfield(P0, 'KP'),     cfg.KP = double(P0.KP); end
        if isfield(P0, 'KI'),     cfg.KI = double(P0.KI); end
        if isfield(P0, 'm'),      cfg.PI_m = double(P0.m); end
        if isfield(P0, 'd'),      cfg.PI_d = double(P0.d); end
    end
    if isfield(pt0, 'DFB') && isstruct(pt0.DFB)
        D0 = pt0.DFB(1);
        if isfield(D0, 'enable'),      cfg.DFB_enable = logical(D0.enable); end
        if isfield(D0, 'gain'),        cfg.DFB_gain = double(D0.gain); end
        if isfield(D0, 'phase_shift'), cfg.DFB_phase_shift = double(D0.phase_shift); end
        if isfield(D0, 'd'),           cfg.DFB_delay = double(D0.d); end
        if isfield(D0, 'm'),           cfg.DFB_m = double(D0.m); end
        if isfield(D0, 'every'),       cfg.DFB_every = double(D0.every); end
    end
    cfg.scan_Df_min = min(Df_scan_kHz);
    cfg.scan_Df_max = max(Df_scan_kHz);
    if length(Df_scan_kHz) > 1, cfg.scan_Df_step = abs(Df_scan_kHz(2) - Df_scan_kHz(1)); end
    cfg.scan_I0_min = min(I0_scan_mA);
    cfg.scan_I0_max = max(I0_scan_mA);
    if length(I0_scan_mA) > 1, cfg.scan_I0_step = abs(I0_scan_mA(2) - I0_scan_mA(1)); end
    sData.config = cfg;
end
save(sumFile, '-struct', 'sData');

% Plot updated stability map
fprintf('Plotting updated stability map...\n');
fig = figure(102); clf;
set(fig, 'Color', 'w', 'Position', [150 150 900 620]);
ax = axes('Parent', fig); hold(ax, 'on');

gr_plot = max(growth_rate, -12.0);
if size(gr_plot, 1) >= 3 && size(gr_plot, 2) >= 3 && all(~isnan(gr_plot(:)))
    gr_pad = [gr_plot(1,1),   gr_plot(1,:),   gr_plot(1,end); ...
              gr_plot(:,1),   gr_plot,        gr_plot(:,end); ...
              gr_plot(end,1), gr_plot(end,:), gr_plot(end,end)];
    k_sm = [0.05 0.10 0.05; 0.10 0.40 0.10; 0.05 0.10 0.05];
    gr_plot = conv2(gr_pad, k_sm, 'valid');
end
growth_rate = gr_plot;

dDf_f = (Df_scan_kHz(2) - Df_scan_kHz(1)) / 2;
dI0_f = (I0_scan_mA(2) - I0_scan_mA(1)) / 2;
Df_fine = linspace(min(Df_scan_kHz) - dDf_f, max(Df_scan_kHz) + dDf_f, 400);
I0_fine = linspace(min(I0_scan_mA) - dI0_f, max(I0_scan_mA) + dI0_f, 300);
[X_orig, Y_orig] = meshgrid(Df_scan_kHz, I0_scan_mA);
[X_fine, Y_fine] = meshgrid(Df_fine, I0_fine);
Xq = max(min(Df_scan_kHz), min(max(Df_scan_kHz), X_fine));
Yq = max(min(I0_scan_mA),  min(max(I0_scan_mA),  Y_fine));
gr_fine = interp2(X_orig, Y_orig, gr_plot, Xq, Yq, 'makima');

imagesc(ax, Df_fine, I0_fine, gr_fine);
set(ax, 'YDir', 'normal');

ylorrd_ctrl = [ ...
    0.00, 0.98, 0.99, 0.97; ...
    0.08, 1.00, 0.98, 0.90; ...
    0.18, 1.00, 0.93, 0.74; ...
    0.32, 1.00, 0.84, 0.56; ...
    0.50, 0.99, 0.70, 0.44; ...
    0.72, 0.96, 0.52, 0.40; ...
    1.00, 0.90, 0.34, 0.32];
u_cmap = linspace(0, 1, 256)';
cmap = max(0, min(1, [interp1(ylorrd_ctrl(:,1), ylorrd_ctrl(:,2), u_cmap, 'pchip'), ...
                      interp1(ylorrd_ctrl(:,1), ylorrd_ctrl(:,3), u_cmap, 'pchip'), ...
                      interp1(ylorrd_ctrl(:,1), ylorrd_ctrl(:,4), u_cmap, 'pchip')]));
colormap(ax, cmap);
valid_gr = growth_rate(~isnan(growth_rate));
if ~isempty(valid_gr) && max(valid_gr) > min(valid_gr)
    clim(ax, [double(min(valid_gr)), double(max(valid_gr))]);
end

cb = colorbar(ax);
cb.Label.String = 'Growth Rate  [s^{-1}]';
cb.Label.FontSize = 11;
cb.Label.FontName = 'Times New Roman';

% Analytical flat-potential curve
I0_line = linspace(min(I0_scan_mA), max(I0_scan_mA), 100) * 1e-3;
try
    Df_opt = arrayfun(@(Ib) detune_HC_calc(Ib, n_hc, C_val, h_val, U0_val, V_mc, R_hc, Q_hc), I0_line) * 1e-3;
    plot(ax, Df_opt, I0_line * 1e3, 'k--', 'LineWidth', 2.0, 'DisplayName', 'Flat-pot. (\xi=1.0)');
catch
end

% Markers
[X_grid, I0_grid] = meshgrid(Df_scan_kHz, I0_scan_mA);
if any(stability_code(:) == 1)
    plot(ax, X_grid(stability_code == 1), I0_grid(stability_code == 1), 'x', 'Color', 'w', 'LineWidth', 3.6, 'MarkerSize', 9, 'HandleVisibility', 'off');
end
plot(ax, X_grid(stability_code == 0), I0_grid(stability_code == 0), 'o', 'MarkerFaceColor', [0.18 0.72 0.22], 'MarkerEdgeColor', [0 0.28 0], 'MarkerSize', 6, 'DisplayName', 'Stable (o)');
plot(ax, X_grid(stability_code == 1), I0_grid(stability_code == 1), 'x', 'Color', 'k', 'LineWidth', 2.0, 'MarkerSize', 8, 'DisplayName', 'PTBL (x)');
plot(ax, X_grid(stability_code == 2), I0_grid(stability_code == 2), 'p', 'MarkerFaceColor', [0.10 0.55 0.98], 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'DisplayName', 'Coupled D-Q (\bigstar)');
plot(ax, X_grid(stability_code == 3), I0_grid(stability_code == 3), '^', 'MarkerFaceColor', [0.58 0.12 0.78], 'MarkerEdgeColor', 'k', 'MarkerSize', 7, 'DisplayName', 'LCBI (\blacktriangle)');
plot(ax, X_grid(stability_code == 4), I0_grid(stability_code == 4), 'd', 'MarkerFaceColor', [0.45 0.45 0.45], 'MarkerEdgeColor', 'k', 'MarkerSize', 6, 'DisplayName', 'Other (\blacklozenge)');

xlabel(ax, 'HHC Detuning \Delta f_{HC}  [kHz]', 'FontSize', 12, 'FontName', 'Times New Roman');
ylabel(ax, 'Beam Current I_0  [mA]', 'FontSize', 12, 'FontName', 'Times New Roman');
title(ax, '2D Stability Map (Full-Bunch Mode Recognition Re-Analyzed)', 'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
legend(ax, 'Location', 'northeast', 'FontSize', 10, 'FontName', 'Times New Roman', 'Box', 'on', 'Color', [1 1 1 0.85]);
grid(ax, 'on'); box(ax, 'on');

out_png = fullfile(folderPath, 'stability_map_full_bunch_reanalyzed.png');
exportgraphics(fig, out_png, 'Resolution', 300);
fprintf('Saved high-res figure to: %s\n', out_png);
fprintf('=================================================================\n');
end
