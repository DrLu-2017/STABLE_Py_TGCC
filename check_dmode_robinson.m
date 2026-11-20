function [is_unstable, Omega_i, thresh_detune, details] = check_dmode_robinson( ...
    I0, R_hc, Q_hc, n_hc, h, f0, U0, V_mc, tau_z, fre_shift_hc, E0, alpha_c, ...
    R_mc, Q_mc, fre_shift_mc, DFB_gain, gr_margin)
% CHECK_DMODE_ROBINSON 快速解析与色散方程联合求解 D-mode / Mode-0 罗宾逊不稳定性
%
% 参考文献:
%   T. He et al., "D-mode Robinson instability in storage rings with a passive
%   superconducting harmonic cavity", Phys. Rev. Accel. Beams 27, 034402 (2024)
%   (arXiv:2311.13205).
%
% 支持 12 个基础参数（向后兼容旧调用）或最多 17 个完整参数（含主腔阻抗、闭环 DFB 等效阻抗与临界缓冲裕度）：
% Inputs:
%   I0           - 束流强度 [A]
%   R_hc, Q_hc   - 谐波腔分路阻抗 [Ohm] 与有载品质因数
%   n_hc         - 谐波次数
%   h            - 主腔谐波数
%   f0           - 回旋频率 [Hz]
%   U0           - 每圈同步辐射能损 [eV]
%   V_mc         - 主腔峰值电压 [V]
%   tau_z        - 纵向辐射阻尼时间 [s]
%   fre_shift_hc - 谐波腔偏调频率 df_rn = f_rn - n*h*f0 [Hz] (束团拉伸时 > 0)
%   E0           - 束流标称能量 [eV]
%   alpha_c      - 动量压实因子
%   R_mc, Q_mc   - (可选) 主腔分路阻抗 [Ohm] 与有载品质因数 (默认 0 = 忽略主腔动态项)
%   fre_shift_mc - (可选) 主腔偏调频率 df_r1 = f_r1 - h*f0 [Hz] (默认 0)
%   DFB_gain     - (可选) 直接射频反馈增益 G_dfb >= 0 (默认 0)
%   gr_margin    - (可选) 预筛选跳过的增长率安全缓冲阈值 [s^-1] (默认 0；预筛选建议设 15~25)
%
% Outputs:
%   is_unstable   - 布尔值，Omega_i > gr_margin 时为 true
%   Omega_i       - 最危险 0 模增长率 [s^-1] (> 0 表示失稳)
%   thresh_detune - D-mode 临界阈值偏调度 df_rn^thresh (Omega_i = 0) [Hz]
%   details       - 结构体，含 f_Dmode_Hz, df1_Hz, Omega_i_Dmode, Omega_i_analytic,
%                   f_Cmode_Hz, Omega_i_Cmode, opt_detune_Hz, is_stretching_safe

    if nargin < 13 || isempty(R_mc),         R_mc = 0;         end
    if nargin < 14 || isempty(Q_mc),         Q_mc = 1;         end
    if nargin < 15 || isempty(fre_shift_mc), fre_shift_mc = 0; end
    if nargin < 16 || isempty(DFB_gain),     DFB_gain = 0;     end
    if nargin < 17 || isempty(gr_margin),    gr_margin = 0;    end

    T0 = 1.0 / f0;
    w0 = 2.0 * pi * f0;

    % 闭环负反馈 DFB 对主腔有效阻抗与品质因数的压低: R1_eff = R1 / (1 + G_dfb / cos(psi1))
    if R_mc > 0 && Q_mc > 0
        f_rf = h * f0;
        psi1 = atan(2.0 * Q_mc * fre_shift_mc / f_rf);
        dfb_fac = 1.0 + max(0.0, DFB_gain) / max(cos(psi1), 0.1);
        R1_eff = R_mc / dfb_fac;
        Q1_eff = Q_mc / dfb_fac;
    else
        R1_eff = 0.0;
        Q1_eff = 1.0;
    end

    % 1. 求解当前工作点的 D-mode 与 C-mode 频率及增长率
    [Om_r_D, Om_i_D, dw1_val, Om_i_he, Om_r_C, Om_i_C] = eval_point_dispersion( ...
        I0, R_hc, Q_hc, n_hc, h, w0, T0, U0, V_mc, tau_z, fre_shift_hc, E0, alpha_c, ...
        R1_eff, Q1_eff, fre_shift_mc);

    % D-mode 专属判据（C-mode 增长率同步保存在 details.Omega_i_Cmode 供诊断）
    Omega_i = Om_i_D;
    is_unstable = (Omega_i > gr_margin);

    % 2. 计算临界阈值偏调度 df_rn^thresh (求解 Omega_i_D(df_hc) = 0, Eq. 26 of He 2024)
    wrn_nom = n_hc * h * w0;
    c_he = (n_hc * h * max(I0, 0) * alpha_c * w0) / (T0 * E0);
    thresh_cubic_Hz = ((2.0 * c_he * R_hc / tau_z)^(1/3)) / (2.0 * pi);

    % 通过对数网格扫描 + 二分搜索精确求解 Omega_i_D(df_hc) = 0 (耗时 < 0.2 ms)
    thresh_detune = find_dmode_threshold_hz( ...
        thresh_cubic_Hz, I0, R_hc, Q_hc, n_hc, h, w0, T0, U0, V_mc, tau_z, E0, alpha_c, ...
        R1_eff, Q1_eff, fre_shift_mc);

    % 3. 计算近最优平坦势束团拉伸偏调度 df_rn^opt (Eq. 29-30 of He 2024)
    k_opt_sq = 1.0 / (n_hc^2) - (1.0 / (n_hc^2 - 1.0)) * (U0 / V_mc)^2;
    if k_opt_sq > 0 && I0 > 0
        Vh_opt = sqrt(k_opt_sq) * V_mc;
        F_opt  = 0.90; % 典型平坦势拉伸束团形状因子幅度 (~0.85-0.95)
        opt_detune_approx_Hz = (F_opt * I0 * wrn_nom * R_hc) / (Vh_opt * Q_hc * 2.0 * pi);
        ratio_hc = (2.0 * I0 * R_hc) / Vh_opt;
        if ratio_hc > 1.0
            opt_detune_flat_Hz = (n_hc * h * f0 / (2.0 * Q_hc)) * sqrt(ratio_hc^2 - 1.0);
        else
            opt_detune_flat_Hz = opt_detune_approx_Hz;
        end
    else
        opt_detune_approx_Hz = 0.0;
        opt_detune_flat_Hz   = 0.0;
    end

    if nargout >= 4
        details = struct();
        details.f_Dmode_Hz          = Om_r_D / (2.0 * pi);
        details.df1_Hz              = dw1_val / (2.0 * pi);
        details.Omega_i_Dmode       = Om_i_D;
        details.Omega_i_analytic    = Om_i_he;
        details.f_Cmode_Hz          = Om_r_C / (2.0 * pi);
        details.Omega_i_Cmode       = Om_i_C;
        details.thresh_cubic_Hz     = thresh_cubic_Hz;
        details.thresh_detune_Hz    = thresh_detune;
        details.opt_detune_Hz       = opt_detune_flat_Hz;
        details.opt_detune_eq30_Hz  = opt_detune_approx_Hz;
        details.is_stretching_safe  = (thresh_detune < opt_detune_flat_Hz);
    end
end

% =========================================================================
% 内部核心求解器：He et al. (PRAB 2024) 闭式解析公式 + 精确复色散方程牛顿修正
% =========================================================================
function [Om_r_D, Om_i_D, delta_w1, Om_i_he, Om_r_C, Om_i_C] = eval_point_dispersion( ...
    I0, Rn, Qn, n_hc, h, w0, T0, U0, Vrf, tau_z, df_rn_Hz, E0, alpha_c, R1, Q1, df_r1_Hz)

    if I0 <= 0 || df_rn_Hz == 0
        w_s0 = sqrt(max(0, h * w0 * alpha_c * sqrt(max(0, Vrf^2 - U0^2)) / (T0 * E0)));
        Om_r_D   = abs(2.0 * pi * df_rn_Hz);
        Om_i_D   = -(n_hc * h * w0) / (2.0 * Qn);
        delta_w1 = 0.0;
        Om_i_he  = Om_i_D;
        Om_r_C   = w_s0;
        Om_i_C   = -1.0 / tau_z;
        return;
    end

    w1_0  = h * w0;
    wn_0  = n_hc * h * w0;
    dw_r1 = 2.0 * pi * df_r1_Hz;
    dw_rn = 2.0 * pi * df_rn_Hz; % 束团拉伸时 df_rn_Hz > 0 => dw_rn > 0
    wr1   = w1_0 + dw_r1;
    wrn   = wn_0 + dw_rn;

    % 计及谐波腔稳态能损修正后的同步相位
    psi_hc   = atan(2.0 * Qn * df_rn_Hz / (wn_0 / (2.0 * pi)));
    U_hc     = 2.0 * I0 * Rn * cos(psi_hc)^2;
    sin_phis = min(0.9999, max(0.0, (U0 + U_hc) / Vrf));
    phi_s    = pi - asin(sin_phis);

    xi    = (h * I0 * alpha_c * w0) / (T0 * E0); % E0 单位为 [eV]，与 Vrf [V] 一致
    ws0_2 = (xi * Vrf * abs(cos(phi_s))) / I0;

    % 谐波腔压平修正后的平衡同步频率平方 (Eq. 5 of He 2024: D-mode 色散与闭式解使用 ws2)
    Zn_stat = cav_impedance(wn_0, Rn, Qn, wrn);
    ws2     = ws0_2 + 2.0 * xi * n_hc * imag(Zn_stat);
    if R1 > 0
        Z1_stat = cav_impedance(w1_0, R1, Q1, wr1);
        ws2_C   = max(1.0, ws2 + 2.0 * xi * imag(Z1_stat));
    else
        ws2_C   = max(1.0, ws2);
    end

    % 1. 主腔相移与边带耦合项 (zeta_1, b_1, k_1)
    if R1 > 0
        psi_minus = atan(2.0 * Q1 * (dw_r1 - dw_rn) / wr1);
        psi_plus  = atan(2.0 * Q1 * (dw_r1 + dw_rn) / wr1);
        zeta1     = -0.5 * xi * R1 * (sin(2.0 * psi_minus) + sin(2.0 * psi_plus));
        b1        = 0.5 * xi * R1 * (cos(2.0 * psi_minus) - cos(2.0 * psi_plus));
        k1        = -(2.0 * xi * R1 * Q1 / wr1) * ...
                    ((cos(psi_minus)^2) * cos(2.0 * psi_minus) - (cos(psi_plus)^2) * cos(2.0 * psi_plus));
    else
        zeta1 = 0.0;
        b1    = 0.0;
        k1    = 0.0;
    end

    % 2. He et al. (PRAB 2024) 一元二次方程闭式解析解
    disc = 1.0;
    if dw_rn > 0
        B1   = (dw_rn / 4.0) + ((zeta1 - ws2) / (4.0 * dw_rn)) - (xi * n_hc * Rn * wrn) / (16.0 * Qn * dw_rn^2);
        B2   = (xi * n_hc * wrn * Rn) / (4.0 * dw_rn * Qn);
        disc = B1^2 - B2;
        if disc >= 0 && B1 > 0
            delta_w1 = B1 - sqrt(disc);
            Om_r_he  = dw_rn - delta_w1;
            bn       = (xi * n_hc * Rn * wrn^2) / (4.0 * Qn^2 * max(delta_w1^2, 1e-18));
            kn       = (xi * n_hc * Rn * wrn)   / (2.0 * Qn   * max(delta_w1^2, 1e-18));
            denom_he = 2.0 * Om_r_he - kn - k1;
            if abs(denom_he) < 1e-9
                denom_he = -1e-9;
            end
            Om_i_he = (bn + b1 - 2.0 * Om_r_he / tau_z) / denom_he;
            Om0_D   = complex(Om_r_he, max(-wrn / Qn, min(1e4, Om_i_he)));
        else
            % 判别式 < 0：进入模式耦合强失稳区
            delta_w1 = real(B1);
            Om_r_he  = max(0.5 * dw_rn, 1.0);
            Om_i_he  = 500.0;
            Om0_D    = complex(0.60 * dw_rn, 0.30 * dw_rn);
        end
    else
        % 负失谐束团缩短区 (Appendix B: dw_rn < 0 -> D-mode 阻尼)
        dw_abs   = abs(dw_rn);
        B1       = (dw_abs / 4.0) - (ws2 / (4.0 * dw_abs)) + (xi * n_hc * Rn * wrn) / (16.0 * Qn * dw_abs^2);
        B2       = (xi * n_hc * wrn * Rn) / (4.0 * dw_abs * Qn);
        delta_w2 = -B1 + sqrt(max(0.0, B1^2 + B2));
        delta_w1 = -delta_w2;
        Om_r_he  = dw_abs + delta_w2;
        bn       = -(xi * n_hc * Rn * wrn^2) / (4.0 * Qn^2 * max(delta_w2^2, 1e-18));
        kn       = -(xi * n_hc * Rn * wrn)   / (2.0 * Qn   * max(delta_w2^2, 1e-18));
        Om_i_he  = (bn + b1 - 2.0 * Om_r_he / tau_z) / (2.0 * Om_r_he - kn - k1);
        Om0_D    = complex(Om_r_he, Om_i_he);
    end

    % 3. 精确复色散方程 g(Omega) = 0 牛顿迭代修正（同时兼容高 Q 超导腔与低 Q 常温铜腔）
    % D-mode 分支（限制在谐波腔极点邻域，防止跳支到 C-mode）
    gam_cav = wrn / (2.0 * Qn);
    dw_abs_ref = abs(dw_rn);
    if gam_cav >= 0.25 * dw_abs_ref
        Om0_D = complex(dw_abs_ref, -gam_cav * 1.005);
    end
    Om_D = Om0_D;
    max_step_D = max([0.35 * abs(delta_w1), 2.0 * gam_cav, 2.0 * pi * 80.0]);
    if disc < 0
        max_step_D = 0.35 * dw_abs_ref;
    end
    for it = 1:20
        [g0, dg] = disp_func_and_deriv(Om_D, w1_0, wn_0, R1, Q1, wr1, Rn, Qn, wrn, xi, n_hc, ws2, tau_z);
        if abs(dg) < 1e-24, break; end
        step = g0 / dg;
        if abs(step) > max_step_D
            step = step * (max_step_D / abs(step));
        end
        Om_next = Om_D - step;
        if abs(real(Om_next) - dw_abs_ref) > 0.65 * dw_abs_ref || real(Om_next) <= 0
            break;
        end
        Om_D = Om_next;
        if abs(step) < 1e-5, break; end
    end

    if abs(real(Om_D) - dw_abs_ref) <= 0.55 * dw_abs_ref && isfinite(Om_D)
        Om_r_D   = abs(real(Om_D));
        Om_i_D   = imag(Om_D);
        delta_w1 = dw_abs_ref - Om_r_D;
    else
        Om_r_D = abs(Om_r_he);
        Om_i_D = Om_i_he;
    end

    % C-mode 分支（同步振荡频率邻域的偶极 Robinson 模式）
    Om_C = complex(sqrt(max(ws0_2, 1e4)), -1.0 / tau_z);
    for it = 1:20
        [g0, dg] = disp_func_and_deriv(Om_C, w1_0, wn_0, R1, Q1, wr1, Rn, Qn, wrn, xi, n_hc, ws2_C, tau_z);
        if abs(dg) < 1e-24, break; end
        step = g0 / dg;
        if abs(step) > 0.4 * abs(Om_C)
            step = step * (0.4 * abs(Om_C) / abs(step));
        end
        Om_C = Om_C - step;
        if abs(step) < 1e-5, break; end
    end
    if isfinite(Om_C) && abs(real(Om_C)) < 0.5 * dw_abs_ref
        Om_r_C = abs(real(Om_C));
        Om_i_C = imag(Om_C);
    else
        Om_r_C = sqrt(max(ws2, 0.0));
        Om_i_C = -1.0 / tau_z;
    end
end

function thresh_hz = find_dmode_threshold_hz( ...
    f_est_hz, I0, Rn, Qn, n_hc, h, w0, T0, U0, Vrf, tau_z, E0, alpha_c, R1, Q1, df_r1_Hz)

    if I0 <= 0 || Rn <= 0 || Qn <= 0
        thresh_hz = 0.0;
        return;
    end

    % 对数网格扫描定位最高 D-mode 失稳过零区间 (100 Hz ~ 600 kHz)
    f_grid = logspace(2.0, 5.78, 40);
    if f_est_hz > 100.0 && f_est_hz < 6e5
        f_grid = sort(unique([f_grid, 0.5*f_est_hz, 0.8*f_est_hz, f_est_hz, 1.25*f_est_hz, 1.6*f_est_hz]));
    end
    gr_grid = zeros(size(f_grid));
    for kg = 1:numel(f_grid)
        [~, gr_grid(kg)] = eval_point_dispersion(I0, Rn, Qn, n_hc, h, w0, T0, U0, Vrf, tau_z, f_grid(kg), E0, alpha_c, R1, Q1, df_r1_Hz);
    end

    pos_idx = find(gr_grid > 0);
    if isempty(pos_idx)
        thresh_hz = 0.0;
        return;
    end
    last_pos = pos_idx(end);
    if last_pos >= numel(f_grid)
        thresh_hz = f_grid(end);
        return;
    end
    f_low  = f_grid(last_pos);
    f_high = f_grid(last_pos + 1);

    % 二分法精确锁定 Omega_i_D = 0 (18 步相对精度 < 1e-5)
    for iter = 1:18
        f_mid = 0.5 * (f_low + f_high);
        [~, gr_mid] = eval_point_dispersion(I0, Rn, Qn, n_hc, h, w0, T0, U0, Vrf, tau_z, f_mid, E0, alpha_c, R1, Q1, df_r1_Hz);
        if gr_mid > 0
            f_low = f_mid;
        else
            f_high = f_mid;
        end
    end
    thresh_hz = 0.5 * (f_low + f_high);
end

function [val, dval] = disp_func_and_deriv(Om, w1_0, wn_0, R1, Q1, wr1, Rn, Qn, wrn, xi, n_hc, ws2, tau_z)
    if R1 > 0
        [Z1p, dZ1p] = cav_impedance(w1_0 + Om, R1, Q1, wr1);
        [Z1m, dZ1m] = cav_impedance(-w1_0 + Om, R1, Q1, wr1);
        dZ1  = Z1p - Z1m;
        ddZ1 = dZ1p - dZ1m;
    else
        dZ1  = 0.0;
        ddZ1 = 0.0;
    end
    [Znp, dZnp] = cav_impedance(wn_0 + Om, Rn, Qn, wrn);
    [Znm, dZnm] = cav_impedance(-wn_0 + Om, Rn, Qn, wrn);
    dZn  = Znp - Znm;
    ddZn = dZnp - dZnm;
    val  = Om^2 + 2i * Om / tau_z - ws2 - 1i * xi * (dZ1 + n_hc * dZn);
    dval = 2.0 * Om + 2i / tau_z - 1i * xi * (ddZ1 + n_hc * ddZn);
end

function [Z, dZ] = cav_impedance(w, R, Q, wr)
    u = Q * (wr / w - w / wr);
    denom = 1.0 + 1i * u;
    if abs(denom) < 1e-14
        denom = 1e-14;
    end
    Z = R / denom;
    if nargout > 1
        du = -Q * (wr / (w * w) + 1.0 / wr);
        dZ = -R * 1i * du / (denom * denom);
    end
end
