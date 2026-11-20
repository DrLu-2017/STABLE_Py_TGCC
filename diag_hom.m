% Diagnostic script for HOM kick
machine_Soleil2;
I0 = 0.5;
fre_shift_hc = 80e3;

fprintf('=== DIAGNOSTIC HOM RUN ===\n');

% Load config from SOLEIL_II_setting.mat
cfg = load('SOLEIL_II_setting.mat');
homData = cfg.config.HOMs;
disp('homData from SOLEIL_II_setting.mat:');
disp(homData);

for test_case = 1:2
    if test_case == 1
        % Table as-is: Col1 = Q, Col2 = R
        fprintf('\n--- TEST CASE 1: Col1=Q, Col2=R ---\n');
        Q_hom   = [cell2mat(homData(:,1))];
        R_hom   = [cell2mat(homData(:,2))];
        fre_hom = [cell2mat(homData(:,3))];
    else
        % Swapped: Col1 = R, Col2 = Q
        fprintf('\n--- TEST CASE 2: Col1=R, Col2=Q (Swapped) ---\n');
        R_hom   = [cell2mat(homData(:,1))];
        Q_hom   = [cell2mat(homData(:,2))];
        fre_hom = [cell2mat(homData(:,3))];
    end
    
    fprintf('Q_hom:   %s\n', mat2str(Q_hom'));
    fprintf('R_hom:   %s Ohm\n', mat2str(R_hom'));
    fprintf('R/Q:     %s Ohm\n', mat2str((R_hom./Q_hom)'));
    fprintf('fre_hom: %s Hz\n', mat2str(fre_hom'));
    
    HALF = machine(C, I0, U0, E0, tau_s, tau_z, sigma_t0, sigma_e0, alpha_c, ...
                   h, V_mc, n_hc, R_hc, Q_hc, fillrate, fre_shift_hc, Q_mc, R_mc, fre_shift_mc);
    HALF.ShortRange_on = 1;
    HALF.Q_hom_m0   = double(Q_hom);
    HALF.R_hom_m0   = double(R_hom);
    HALF.fre_hom_m0 = double(fre_hom);
    HALF.wrf_hom_m0       = 2*pi*HALF.fre_hom_m0;
    HALF.angle_hom_m0     = HALF.wrf_hom_m0 * HALF.Tb;
    HALF.VbImagFactor_hom_m0 = 1./sqrt(4*HALF.Q_hom_m0.^2-1);
    HALF.rot_coef_hom_m0     = sqrt(1-1./(2*HALF.Q_hom_m0).^2);
    HALF.Q_hom_length = length(HALF.Q_hom_m0);
    
    Bun_num = sum(pattern==1);
    Par_num = 1000;
    
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
    
    rng(42);
    q1 = TruncatedGaussian(1, [-3,3], [Par_num, 1]);
    p1 = TruncatedGaussian(1, [-3,3], [Par_num, 1]);
    Q = gpuArray(single(repmat(q1, 1, Bun_num)));
    P = gpuArray(single(repmat(p1, 1, Bun_num)));
    index_add = gpuArray(single((1:Bun_num) - 1));
    
    V_hc_load_0 = HALF.V_hc_load_0;
    V_mc_load_0 = HALF.V_mc_load_0;
    HALF.Vg_mc_track_0 = HALF.Vg_mc_init;
    PI = []; DFB = [];
    
    for i = 1:10
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
        fprintf('  Turn %2d: binnum = %g, max|P| = %.2e, max|Q| = %.2e\n', ...
            i, binnum, gather(max(max(abs(P)))), gather(max(max(abs(Q)))));
        if binnum > 2000 || isnan(binnum) || isinf(binnum) || binnum < 1
            fprintf('  -> BEAM LOST on turn %d!\n', i);
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
    end
end
exit;
