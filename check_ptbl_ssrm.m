function is_ptbl_unstable = check_ptbl_ssrm(I0, R_hc, Q_hc, n_hc, h, f0, U0, V_mc, fre_shift_hc, fillrate, F_hc)
% CHECK_PTBL_SSRM 快速解析预筛选深度过拉伸双势阱 PTBL 失稳区
%
% 说明：仅对远离平坦势临界区 (xi_eff > 1.15 且高流强) 的确定性深度双势阱过拉伸失稳点返回 true，
%       保留临界区 (0.85 <= xi_eff <= 1.15) 交由完整 STABLE 多粒子 GPU 追踪精确评估朗道阻尼与闭环效应。

    if nargin < 10 || isempty(fillrate)
        fillrate = 1.0;
    end
    if nargin < 11 || isempty(F_hc)
        F_hc = 1.0;
    end

    if I0 <= 0 || fre_shift_hc <= 0
        is_ptbl_unstable = false;
        return;
    end

    f_hc   = n_hc * h * f0;
    psi_hc = atan(2.0 * Q_hc * fre_shift_hc / f_hc);
    I_eff  = I0 * min(max(fillrate, 0.1), 1.0);
    V_hc   = 2.0 * F_hc * I_eff * R_hc * cos(psi_hc);

    k_fp_sq = 1.0 / (n_hc^2) - (1.0 / (n_hc^2 - 1.0)) * (U0 / V_mc)^2;
    if k_fp_sq <= 0
        is_ptbl_unstable = false;
        return;
    end
    k_fp = sqrt(k_fp_sq);
    xi   = V_hc / (k_fp * V_mc);

    % 保守预筛阈值：仅跳过深度双势阱过拉伸区 (xi_eff > 1.15 且 I0 >= 365 mA)
    is_ptbl_unstable = (xi > 1.15) && (I0 >= 0.365);
end
