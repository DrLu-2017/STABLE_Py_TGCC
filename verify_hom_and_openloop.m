% Verification script for open-loop cavity tracking, HOM tracking, and classification
fprintf('========================================================\n');
fprintf('  VERIFICATION: Open-loop MC generator drive & HOM test \n');
fprintf('========================================================\n');

machine_Soleil2;
I0 = 0.3; % 300 mA
fre_shift_hc = 80e3;

HALF = machine(C, I0, U0, E0, tau_s, tau_z, sigma_t0, sigma_e0, alpha_c, ...
               h, V_mc, n_hc, R_hc, Q_hc, fillrate, fre_shift_hc, Q_mc, R_mc, fre_shift_mc);
HALF.ShortRange_on = 1;

% HOMs setup (Soleil 2 3-HOM setup)
Q_hom_m0   = [830; 921; 209];
R_hom_m0   = [6500; 7300; 10100];
fre_hom_m0 = [1.6971e9; 1.18e9; 959.1e6];

HALF.Q_hom_m0   = Q_hom_m0;
HALF.R_hom_m0   = R_hom_m0;
HALF.fre_hom_m0 = fre_hom_m0;
HALF.wrf_hom_m0       = 2*pi*HALF.fre_hom_m0;
HALF.angle_hom_m0     = HALF.wrf_hom_m0 * HALF.Tb;
HALF.VbImagFactor_hom_m0 = 1./sqrt(4*HALF.Q_hom_m0.^2-1);
HALF.rot_coef_hom_m0     = sqrt(1-1./(2*HALF.Q_hom_m0).^2);
HALF.Q_hom_length = length(HALF.Q_hom_m0);

Bun_num = sum(pattern==1);
Par_num = 1000;
Track_num = 500;
Recor_step = 10;
Recor_num  = floor(Track_num / Recor_step);

charge = ones(1, h) .* pattern;
charge = charge / sum(charge) * Bun_num;
HALF.qc    = charge .* pattern * HALF.qc / Par_num;
HALF.Vb_hc = HALF.qc * HALF.wr_hc * R_hc / Q_hc;
HALF.Vb_mc = HALF.qc * HALF.wr_mc * R_mc / Q_mc;

Dq = 0.05;
rot_decay_coef_hc = 1i*HALF.rot_coef_hc - 1/(2*Q_hc);
TbAng_coef_hc     = exp(rot_decay_coef_hc * HALF.angle_hc);
exp_ang_coef_hc   = -rot_decay_coef_hc * HALF.wr_hc * sigma_t0;

rot_decay_coef_mc = 1i*HALF.rot_coef_mc - 1/(2*Q_mc);
TbAng_coef_mc     = exp(rot_decay_coef_mc * HALF.angle_mc);
exp_ang_coef_mc   = -rot_decay_coef_mc * HALF.wr_mc * sigma_t0;

HALF.Vb_hom = HALF.qc .* (HALF.wrf_hom_m0 .* HALF.R_hom_m0 ./ HALF.Q_hom_m0);
rot_decay_coef_hom = 1i * HALF.rot_coef_hom_m0 - 1 ./ (2 * HALF.Q_hom_m0);
TbAng_coef_hom     = exp(rot_decay_coef_hom .* HALF.angle_hom_m0);
exp_ang_coef_hom   = -rot_decay_coef_hom .* HALF.wrf_hom_m0 * sigma_t0;
V_hom_load_0       = zeros(HALF.Q_hom_length, 1) + 0i;

wake_kick_coef = HALF.qc * HALF.kick_coef;
tau_q = (0:Dq:200)' * sigma_t0;
Wake_inter_hc = -HALF.wr_hc * R_hc/Q_hc * exp(-tau_q*HALF.wr_hc/(2*Q_hc)) .* ...
    (cos(tau_q*HALF.wr_hc*HALF.rot_coef_hc) - HALF.VbImagFactor_hc * sin(tau_q*HALF.wr_hc*HALF.rot_coef_hc));
Wake_inter_hc(1) = Wake_inter_hc(1)/2;
Wake_inter_mc = -HALF.wr_mc * R_mc/Q_mc * exp(-tau_q*HALF.wr_mc/(2*Q_mc)) .* ...
    (cos(tau_q*HALF.wr_mc*HALF.rot_coef_mc) - HALF.VbImagFactor_mc * sin(tau_q*HALF.wr_mc*HALF.rot_coef_mc));
Wake_inter_mc(1) = Wake_inter_mc(1)/2;
Wake_inter = Wake_inter_hc + Wake_inter_mc;
for jj = 1:length(HALF.Q_hom_m0)
    Wake_inter_hom = -HALF.wrf_hom_m0(jj)*HALF.R_hom_m0(jj)/HALF.Q_hom_m0(jj) * ...
        exp(-tau_q*HALF.wrf_hom_m0(jj)/2/HALF.Q_hom_m0(jj)) .* ...
        (cos(tau_q*HALF.wrf_hom_m0(jj)*HALF.rot_coef_hom_m0(jj)) - ...
         HALF.VbImagFactor_hom_m0(jj)*sin(tau_q*HALF.wrf_hom_m0(jj)*HALF.rot_coef_hom_m0(jj)));
    Wake_inter_hom(1) = Wake_inter_hom(1)/2;
    Wake_inter = Wake_inter + Wake_inter_hom;
end
Wake_inter = gpuArray(single(Wake_inter));

% Test Case: PI feedback DISABLED (open loop), DFB DISABLED
PI = struct();
PI.m         = 5120;
PI.d         = 800;
PI.piIndex   = 0;
PI.Integral  = 0;
PI.RL        = HALF.R_mc;
PI.Ts        = PI.m * HALF.Tb;
PI.KP        = 0;
PI.KI        = 0;
PI.Ig0       = HALF.Ig_mc_0;
PI.DIg       = 0;
PI.Ig        = PI.Ig0;
PI.Vg0       = HALF.Vg_mc_0;
PI.DI        = 0;
PI.Ig_track  = [];
PI.Ig_FBid   = [];
PI.RoverQ    = HALF.R_mc / HALF.Q_mc;
PI.Vrf_mc_track = zeros(1, PI.m);
PI.enable    = false; % Open loop!
PI.Enabled   = false;
DFB = [];

rng(42);
q1 = TruncatedGaussian(1, [-3,3], [Par_num, 1]);
p1 = TruncatedGaussian(1, [-3,3], [Par_num, 1]);
Q = gpuArray(single(repmat(q1, 1, Bun_num)));
P = gpuArray(single(repmat(p1, 1, Bun_num)));
index_add = gpuArray(single((1:Bun_num) - 1));

V_hc_load_0 = HALF.V_hc_load_0;
V_mc_load_0 = HALF.V_mc_load_0;
HALF.Vg_mc_track_0 = HALF.Vg_mc_init;
beam_lost = false;

record_Q_mean = zeros(Recor_num, Bun_num);
record_Q_std  = zeros(Recor_num, Bun_num);
record_P_mean = zeros(Recor_num, Bun_num);
record_P_std  = zeros(Recor_num, Bun_num);
record_th = 0;

for i = 1:Track_num
    Q = Q + P * HALF.drift_coef;
    Q_min = min(Q);
    Q_new = round((Q - Q_min) * (1/Dq));

    for jj = 1:HALF.Q_hom_length
        exp_angle = exp(exp_ang_coef_hom(jj) * Q);
        exp_angle_sum = gather(sum(exp_angle));
        V_load_cpu = double(exp_angle_sum .* HALF.Vb_hom(jj, pattern==1));
        [V_load, V_hom_load_0(jj)] = VoltageLoadCalc_matlab(V_hom_load_0(jj), V_load_cpu, TbAng_coef_hom(jj), pattern);
        V_load_cpu = V_load(pattern==1) * HALF.kick_coef;
        V_load_gpu = gpuArray(single(V_load_cpu));
        V_hom_load_kick = V_load_gpu ./ exp_angle;
        P = P - real(V_hom_load_kick) + imag(V_hom_load_kick) * HALF.VbImagFactor_hom_m0(jj);
    end

    exp_angle = exp(exp_ang_coef_hc * Q);
    exp_angle_sum = gather(sum(exp_angle));
    V_load_cpu = double(exp_angle_sum .* HALF.Vb_hc(pattern==1));
    [V_load, V_hc_load_0] = VoltageLoadCalc_matlab(V_hc_load_0, V_load_cpu, TbAng_coef_hc, pattern);
    V_load_cpu = V_load(pattern==1) * HALF.kick_coef;
    V_hc_load_kick = gpuArray(single(V_load_cpu)) ./ exp_angle;

    exp_angle = exp(exp_ang_coef_mc * Q);
    exp_angle_sum = gather(sum(exp_angle));
    V_load_cpu = double(exp_angle_sum .* HALF.Vb_mc(pattern==1));
    [Vc_mc, ~, HALF.Vg_mc_track_0, V_load, HALF.V_mc_load_0, PI, DFB] = ...
        PI_Control(PI, HALF.Vrf_ideal, HALF.Vg_mc_track_0, HALF.V_mc_load_0, V_load_cpu, TbAng_coef_mc, pattern, DFB);
    V_mc_kick = gpuArray(single(Vc_mc(pattern==1) * HALF.kick_coef)) ./ exp_angle;

    binnum = gather(max(max(Q_new))) + 1;
    if binnum > 2000 || isnan(binnum) || isinf(binnum) || binnum < 1
        beam_lost = true;
        fprintf('  FAILED: Beam lost at turn %d! binnum = %g\n', i, binnum);
        break;
    end

    bin_num_q = sum(BinNumCalZ(binnum, Q_new));
    bin_num_q = reshape(bin_num_q, binnum, Bun_num);
    kick_conv = conv2(bin_num_q, Wake_inter(1:binnum));
    kick_conv = kick_conv(1:binnum,:) .* wake_kick_coef(pattern==1) * min(i/5000, 1);
    Q_new = Q_new + (1 + index_add * binnum);
    wake_kick = kick_conv(Q_new);

    rad_quan_kick = -HALF.radampcoef * P + HALF.quanexcoef * gpuArray.randn(Par_num, Bun_num, 'single');
    P = P + rad_quan_kick - HALF.ploss - real(V_hc_load_kick) + imag(V_hc_load_kick)*HALF.VbImagFactor_hc ...
        - real(V_mc_kick) + imag(V_mc_kick)*HALF.VbImagFactor_mc + wake_kick;

    if mod(i, Recor_step) == 0
        record_th = record_th + 1;
        record_Q_mean(record_th,:) = gather(mean(Q));
        record_Q_std(record_th,:)  = gather(std(Q));
        record_P_mean(record_th,:) = gather(mean(P));
        record_P_std(record_th,:)  = gather(std(P));
    end
end

if ~beam_lost
    fprintf('  SUCCESS: Completed %d turns without beam loss!\n', Track_num);
    fprintf('  Final mean|P| = %.4f, mean|Q| = %.4f, Vg_mc_track_0 = %.2f MV\n', ...
        mean(abs(record_P_mean(end,:))), mean(abs(record_Q_mean(end,:))), abs(HALF.Vg_mc_track_0)*1e-6);
    
    % Test classification logic:
    win_analysis = round(0.5*Recor_num)+1 : Recor_num;
    centroid_E = record_P_mean(win_analysis, 1)';
    centroid_Q = record_Q_mean(win_analysis, 1)';
    sigma_BL   = record_Q_std(win_analysis, :) * sigma_t0;
    
    xx = log(abs(centroid_E) + 1e-30);
    pass_turns = win_analysis * Recor_step;
    a_fit = polyfit(pass_turns, xx, 1);
    gr_si = a_fit(1) / HALF.T0;
    
    sig_fft = centroid_Q - mean(centroid_Q);
    N_fft = length(sig_fft);
    fft_amp = abs(fft(sig_fft));
    fft_amp = fft_amp(1:floor(N_fft/2));
    freq_axis = (0:floor(N_fft/2)-1) / N_fft / (Recor_step * HALF.T0);
    [~, peak_idx] = max(fft_amp(2:end));
    f_peak_Hz = freq_axis(peak_idx + 1);
    
    BL_mean_all = mean(sigma_BL, 2);
    BL_osc_amp = std(BL_mean_all);
    BL_osc_rel = BL_osc_amp / mean(BL_mean_all);
    
    Q_mean_last = mean(record_Q_mean(end-min(99, size(record_Q_mean,1)-1):end, :), 1);
    std_bb = std(Q_mean_last) * sigma_t0 * 1e12;
    
    win_steady = max(1, round(0.8 * Recor_num)) : Recor_num;
    bunch_jitters = std(record_Q_mean(win_steady, :), 0, 1) * sigma_t0 * 1e12;
    centroid_jitter_ps = mean(bunch_jitters);
    max_bunch_jitter_ps = max(bunch_jitters);
    
    damping_rate = 1 / HALF.tau_z;
    noise_floor_ps = (sigma_t0 * 1e12) / sqrt(max(Par_num, 100));
    jitter_thresh_mean = max(1.2, 3.0 * noise_floor_ps);
    jitter_thresh_max  = max(2.0, 5.0 * noise_floor_ps);
    
    is_unstable = (gr_si > 0.3 * damping_rate) || ...
                  (centroid_jitter_ps > jitter_thresh_mean) || ...
                  (max_bunch_jitter_ps > jitter_thresh_max && I0 >= 0.35) || ...
                  (BL_osc_rel > 0.05) || ...
                  (std_bb > 0.6);
                  
    psi_hc_calc = atan(2 * Q_hc * fre_shift_hc / HALF.fre_hc);
    V_hc_calc   = 2 * I0 * R_hc * cos(psi_hc_calc);
    k_fp_calc   = sqrt(1/n_hc^2 - 1/(n_hc^2-1) * (U0/V_mc)^2);
    xi_calc     = V_hc_calc / (k_fp_calc * V_mc);
    
    code = 0;
    if is_unstable
        if (xi_calc >= 1.05 || (std_bb > 0.6 && mean(BL_mean_all)*sigma_t0*1e12 > 75)) && I0 >= 0.375
            code = 1;
        elseif ~isempty(Q_hom_m0) && xi_calc < 0.98 && ...
               (f_peak_Hz > 2 * (HALF.f_rf*2.05e-3/416) || gr_si > 0.5 * damping_rate || centroid_jitter_ps > jitter_thresh_mean)
            code = 3;
        elseif BL_osc_rel > 0.05 || centroid_jitter_ps > 1.0 || (xi_calc >= 0.85 && xi_calc <= 1.05)
            code = 2;
        else
            code = 4;
        end
    end
    fprintf('  Stability classification code = %d (Expected: 0 for Stable)\n', code);
    fprintf('  gr_si = %.2f s^-1, centroid_jitter = %.3f ps, noise_floor = %.3f ps\n', ...
        gr_si, centroid_jitter_ps, noise_floor_ps);
    fprintf('========================================================\n');
end
