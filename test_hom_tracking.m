% Test HOM tracking and stability classification
machine_Soleil2;
I0 = 0.5;
fre_shift_hc = 80e3;

% Case 1: as parsed currently (Col 1 = Q, Col 2 = R)
% homData = {830, 5500, 1.6971e9; 521, 7300, 1.18e9; 209, 10100, 959.1e6};
% Case 2: if Col 1 = R, Col 2 = Q
% homData = {5500, 830, 1.6971e9; 7300, 521, 1.18e9; 10100, 209, 959.1e6};

for test_case = 1:2
    if test_case == 1
        fprintf('\n=== Test Case 1: Col1 = Q (830), Col2 = R (5500) ===\n');
        Q_hom_m0   = [830; 521; 209];
        R_hom_m0   = [5500; 7300; 10100];
        fre_hom_m0 = [1.6971e9; 1.18e9; 959.1e6];
    else
        fprintf('\n=== Test Case 2: Col1 = R (830), Col2 = Q (5500) ===\n');
        Q_hom_m0   = [5500; 7300; 10100];
        R_hom_m0   = [830; 521; 209];
        fre_hom_m0 = [1.6971e9; 1.18e9; 959.1e6];
    end
    
    HALF = machine(C, I0, U0, E0, tau_s, tau_z, sigma_t0, sigma_e0, alpha_c, ...
                   h, V_mc, n_hc, R_hc, Q_hc, fillrate, fre_shift_hc, Q_mc, R_mc, fre_shift_mc);
    HALF.ShortRange_on = 1;
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
    Track_num = 200;
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

    q1 = randn(Par_num, 1);
    p1 = randn(Par_num, 1);
    Q = gpuArray(single(repmat(q1, 1, Bun_num)));
    P = gpuArray(single(repmat(p1, 1, Bun_num)));
    index_add = gpuArray(single((1:Bun_num) - 1));

    V_hc_load_0 = HALF.V_hc_load_0;
    V_mc_load_0 = HALF.V_mc_load_0;
    HALF.Vg_mc_track_0 = HALF.Vg_mc_init;
    PI = []; DFB = [];
    beam_lost = false;
    record_th = 0;
    record_Q_mean = zeros(Recor_num, Bun_num);
    record_Q_std  = zeros(Recor_num, Bun_num);
    record_P_mean = zeros(Recor_num, Bun_num);
    record_P_std  = zeros(Recor_num, Bun_num);

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
            fprintf('  Turn %d: BEAM LOST! binnum = %g\n', i, binnum);
            break;
        end
        try
            bin_num_q = sum(BinNumCalZ(binnum, Q_new));
        catch ME
            beam_lost = true;
            fprintf('  Turn %d: CUDA EXCEPTION! %s\n', i, ME.message);
            break;
        end
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
        fprintf('  Tracking finished %d turns. Max binnum = %d\n', Track_num, binnum);
        win_analysis = round(0.5*Recor_num)+1 : Recor_num;
        centroid_E = record_P_mean(win_analysis, 1)';
        centroid_Q = record_Q_mean(win_analysis, 1)';
        sigma_BL   = record_Q_std(win_analysis, :) * sigma_t0;
        xx = log(abs(centroid_E) + 1e-30);
        pass_turns = win_analysis * Recor_step;
        a_fit = polyfit(pass_turns, xx, 1);
        gr_si = a_fit(1) / HALF.T0;
        BL_mean_all = mean(sigma_BL, 2);
        BL_osc_amp = std(BL_mean_all);
        BL_osc_rel = BL_osc_amp / mean(BL_mean_all);
        bunch_jitters = std(record_Q_mean(win_analysis, :), 0, 1) * sigma_t0 * 1e12;
        centroid_jitter_ps = mean(bunch_jitters);
        max_bunch_jitter_ps = max(bunch_jitters);
        std_bb = std(mean(record_Q_mean(win_analysis, :), 1)) * sigma_t0 * 1e12;
        
        damping_rate = 1 / HALF.tau_z;
        fprintf('  gr_si = %.2e (0.3*damping = %.2e)\n', gr_si, 0.3*damping_rate);
        fprintf('  centroid_jitter_ps = %.3f ps (threshold = 1.2)\n', centroid_jitter_ps);
        fprintf('  max_bunch_jitter_ps = %.3f ps (threshold = 0.30)\n', max_bunch_jitter_ps);
        fprintf('  BL_osc_rel = %.4f (threshold = 0.05)\n', BL_osc_rel);
        fprintf('  std_bb = %.3f ps (threshold = 0.6)\n', std_bb);

        is_unstable = (gr_si > 0.3 * damping_rate) || ...
                      (centroid_jitter_ps > 1.2) || ...
                      (max_bunch_jitter_ps > 0.30 && I0 >= 0.35) || ...
                      (BL_osc_rel > 0.05) || ...
                      (std_bb > 0.6);
        fprintf('  is_unstable = %d\n', is_unstable);

        psi_hc_calc = atan(2 * Q_hc * fre_shift_hc / HALF.fre_hc);
        V_hc_calc   = 2 * I0 * R_hc * cos(psi_hc_calc);
        k_fp_calc   = sqrt(1/n_hc^2 - 1/(n_hc^2-1) * (U0/V_mc)^2);
        xi_calc     = V_hc_calc / (k_fp_calc * V_mc);
        fprintf('  xi_calc = %.3f\n', xi_calc);

        code = 0;
        if is_unstable
            if (xi_calc >= 1.05 || (std_bb > 0.6 && mean(BL_mean_all)*sigma_t0*1e12 > 75)) && I0 >= 0.375
                code = 1;
            elseif xi_calc < 0.98 && (fre_shift_hc*1e-3 >= 110 || (max_bunch_jitter_ps > 0.28 && fre_shift_hc*1e-3 >= 90))
                code = 3;
            elseif BL_osc_rel > 0.05 || centroid_jitter_ps > 1.0 || (xi_calc >= 0.85 && xi_calc <= 1.05)
                code = 2;
            else
                code = 4;
            end
        end
        fprintf('  Classification code = %d\n', code);
    end
end
