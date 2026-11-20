function out = calc_bunch_form_factor_eq(mode, I0, V_mc, U0, E0, C, h, alpha_c, sigma_e0, n_hc, R_hc, Q_hc, param_val, use_form_factor)
% CALC_BUNCH_FORM_FACTOR_EQ Self-consistent Fokker-Planck potential-well equilibrium
% solver for passive harmonic cavity bunch form factor F_n(sigma_tau) & flat-potential detuning.
%
% Modes:
%   1) mode = 'eval':
%      param_val = fre_shift_hc [Hz]
%      Computes self-consistent equilibrium F_hc, F_mc, sigma_t_eq, fais_mc_whc,
%      V_hc_load_0, and effective flat-potential parameter xi_eff at given detuning.
%
%   2) mode = 'detune_for_xi':
%      param_val = xi_target (e.g. 1.0, 0.98, 1.05)
%      Computes self-consistent HHC detuning fre_shift_hc [Hz] and tuning angle psi_hc [deg]
%      such that 2 * F_hc(sigma_tau) * I0 * R_hc * cos(psi_hc) / (k_fp * V_mc) = xi_target.

if nargin < 14 || isempty(use_form_factor)
    use_form_factor = true;
end

cspeed = 299792458.0;
T0     = C / cspeed;
f0     = 1.0 / T0;
f_rf   = h * f0;
w_rf   = 2.0 * pi * f_rf;
f_hc   = n_hc * f_rf;
wre_hc = 2.0 * pi * f_hc;

k_fp_sq = 1.0 / (n_hc^2) - (1.0 / max(n_hc^2 - 1.0, 1.0)) * (U0 / max(V_mc, 1.0))^2;
k_fp    = sqrt(max(k_fp_sq, 1e-8));

Qs0      = sqrt(max(0.0, alpha_c / (2.0 * pi) * h * sqrt(max(0.0, (V_mc / E0)^2 - (U0 / E0)^2))));
sigma_t0 = alpha_c * sigma_e0 / max(2.0 * pi * Qs0, 1e-12) * T0;

out = struct();
out.sigma_t0 = sigma_t0;
out.k_fp     = k_fp;

if strcmpi(mode, 'detune_for_xi')
    xi_target = param_val;
    if I0 <= 1e-9 || R_hc <= 0 || Q_hc <= 0 || k_fp <= 0
        out.fre_shift_hc = 0.0;
        out.psi_hc_deg   = 0.0;
        out.F_hc         = 1.0;
        out.F_mc         = 1.0;
        out.sigma_t_eq   = sigma_t0;
        out.xi_eff       = 0.0;
        return;
    end

    if ~use_form_factor
        cos_psi = (xi_target * k_fp * V_mc) / (2.0 * I0 * R_hc);
        if cos_psi >= 1.0
            psi = 0.0;
            df_hc = 0.0;
        else
            psi = acos(max(-1.0, cos_psi));
            df_hc = tan(psi) * f_hc / (2.0 * Q_hc);
        end
        out.fre_shift_hc = df_hc;
        out.psi_hc_deg   = psi * 180.0 / pi;
        out.F_hc         = 1.0;
        out.F_mc         = 1.0;
        out.sigma_t_eq   = sigma_t0;
        out.xi_eff       = (2.0 * I0 * R_hc * cos(psi)) / (k_fp * V_mc);
        return;
    end

    % Self-consistent Fokker-Planck potential-well iteration along constant xi_eff = xi_target
    tau_max = min(0.45 / f_rf, max(12.0 * sigma_t0, 320e-12));
    tau     = linspace(-tau_max, tau_max, 401);
    dtau    = tau(2) - tau(1);

    Fn  = 0.92;
    F1  = 0.995;
    th1 = 0.0;
    thn = 0.0;
    psi = 0.0;
    sig_t = sigma_t0;

    for it = 1:25
        cos_psi = (xi_target * k_fp * V_mc) / max(2.0 * Fn * I0 * R_hc, 1e-9);
        if cos_psi >= 0.999999
            psi = 0.0;
        else
            psi = acos(max(0.0, cos_psi));
        end

        arg_s = (U0 + 2.0 * (Fn^2) * I0 * R_hc * (cos(psi)^2)) / max(F1 * V_mc, 1.0);
        arg_s = max(-0.999, min(0.999, arg_s));
        phis  = pi - asin(arg_s) - th1;

        Fn_c   = Fn * exp(-1i * thn);
        E_gain = V_mc * sin(phis + w_rf * tau) ...
               - 2.0 * I0 * R_hc * cos(psi) * real(exp(1i * psi) * Fn_c * exp(1i * n_hc * w_rf * tau)) ...
               - U0;

        U_pot = -[0.0, cumsum(0.5 * (E_gain(1:end-1) + E_gain(2:end))) * dtau] / (T0 * E0);
        U_pot = U_pot - min(U_pot);

        exponent = min(100.0, max(0.0, U_pot / max(alpha_c * sigma_e0^2, 1e-18)));
        rho = exp(-exponent);
        rho = rho / max(sum(rho) * dtau, 1e-30);

        t_mean = sum(tau .* rho) * dtau;
        tau_c  = tau - t_mean;
        sig_t  = sqrt(sum((tau_c.^2) .* rho) * dtau);

        F1_new = sum(rho .* exp(-1i * w_rf * tau_c)) * dtau;
        Fn_new = sum(rho .* exp(-1i * n_hc * w_rf * tau_c)) * dtau;

        err = abs(abs(Fn_new) - Fn);
        Fn  = 0.5 * abs(Fn_new) + 0.5 * Fn;
        F1  = 0.5 * abs(F1_new) + 0.5 * F1;
        th1 = -angle(F1_new);
        thn = -angle(Fn_new);

        if err < 1e-6 && it >= 4
            break;
        end
    end

    df_hc = tan(psi) * f_hc / (2.0 * Q_hc);
    out.fre_shift_hc = df_hc;
    out.psi_hc_deg   = psi * 180.0 / pi;
    out.F_hc         = Fn;
    out.F_mc         = F1;
    out.sigma_t_eq   = sig_t;
    out.xi_eff       = (2.0 * Fn * I0 * R_hc * cos(psi)) / (k_fp * V_mc);
    return;
end

% Mode = 'eval' at specified fre_shift_hc [Hz]
fre_shift_hc = param_val;
wr_hc = wre_hc + 2.0 * pi * fre_shift_hc;
det_angle_hc = pi - atan(Q_hc * (wr_hc / wre_hc - wre_hc / wr_hc));
psi_hc = pi - det_angle_hc; % in (-pi/2, pi/2)

if ~use_form_factor || I0 <= 1e-9 || R_hc <= 0 || Q_hc <= 0
    arg_s = max(-0.9999, min(0.9999, (U0 + 2.0 * I0 * R_hc * cos(det_angle_hc)^2) / max(V_mc, 1.0)));
    out.F_hc         = 1.0;
    out.F_mc         = 1.0;
    out.theta_hc     = 0.0;
    out.theta_mc     = 0.0;
    out.sigma_t_eq   = sigma_t0;
    out.det_angle_hc = det_angle_hc;
    out.psi_hc_deg   = psi_hc * 180.0 / pi;
    out.U_hc_real    = 2.0 * I0 * R_hc * cos(psi_hc)^2;
    out.fais_mc_whc  = pi - asin(arg_s);
    out.V_hc_load_0  = -2.0 * I0 * R_hc * cos(det_angle_hc) * exp(1i * (pi - det_angle_hc));
    out.xi_eff       = (2.0 * I0 * R_hc * abs(cos(psi_hc))) / max(k_fp * V_mc, 1.0);
    return;
end

tau_max = min(0.45 / f_rf, max(12.0 * sigma_t0, 320e-12));
tau     = linspace(-tau_max, tau_max, 401);
dtau    = tau(2) - tau(1);

Fn_c = 0.90;
F1_c = 0.994;
sig_t = sigma_t0;
converged = false;
phis_cent = pi - asin(min(0.999, U0 / max(V_mc, 1.0)));

for it = 1:35
    arg_s = (U0 + 2.0 * (abs(Fn_c)^2) * I0 * R_hc * (cos(psi_hc)^2)) / max(abs(F1_c) * V_mc, 1.0);
    arg_s = max(-0.999, min(0.999, arg_s));
    phis  = pi - asin(arg_s) + angle(F1_c);

    E_gain = V_mc * sin(phis + w_rf * tau) ...
           - 2.0 * I0 * R_hc * cos(psi_hc) * real(exp(1i * psi_hc) * Fn_c * exp(1i * n_hc * w_rf * tau)) ...
           - U0;

    U_pot = -[0.0, cumsum(0.5 * (E_gain(1:end-1) + E_gain(2:end))) * dtau] / (T0 * E0);
    U_pot = U_pot - min(U_pot);

    exponent = min(100.0, max(0.0, U_pot / max(alpha_c * sigma_e0^2, 1e-18)));
    rho = exp(-exponent);
    rho = rho / max(sum(rho) * dtau, 1e-30);

    t_mean = sum(tau .* rho) * dtau;
    tau_c  = tau - t_mean;
    sig_t  = sqrt(sum((tau_c.^2) .* rho) * dtau);

    F1_new = sum(rho .* exp(-1i * w_rf * tau_c)) * dtau;
    Fn_new = sum(rho .* exp(-1i * n_hc * w_rf * tau_c)) * dtau;

    err  = abs(Fn_new - Fn_c);
    Fn_c = 0.35 * Fn_new + 0.65 * Fn_c;
    F1_c = 0.35 * F1_new + 0.65 * F1_c;
    phis_cent = phis + w_rf * t_mean;

    if err < 1e-6 && it >= 5
        converged = true;
        break;
    end
end

% Guard against over-stretched double-well asymmetric bucket jump (PTBL regime, xi_eff > 1.05):
% If stationary single-bunch equilibrium bifurcated (phis_cent > pi - 0.03 or non-convergent sawtooth),
% evaluate smooth flat-potential envelope at xi_target = 1.04 so initial synchronous phase stays physical.
xi_raw = (2.0 * abs(Fn_c) * I0 * R_hc * abs(cos(psi_hc))) / max(k_fp * V_mc, 1.0);
if ~converged || phis_cent > (pi - 0.03) || phis_cent < (0.5 * pi) || (xi_raw > 1.06 && sig_t < 1.8 * sigma_t0)
    fp_ref = calc_bunch_form_factor_eq('detune_for_xi', I0, V_mc, U0, E0, C, h, alpha_c, sigma_e0, n_hc, R_hc, Q_hc, 1.04, true);
    Fn_mag = fp_ref.F_hc;
    F1_mag = fp_ref.F_mc;
    th_hc  = 0.0;
    sig_t  = fp_ref.sigma_t_eq;
    arg_s  = max(-0.999, min(0.999, (U0 + 2.0 * (Fn_mag^2) * I0 * R_hc * cos(psi_hc)^2) / max(F1_mag * V_mc, 1.0)));
    phis_cent = pi - asin(arg_s);
else
    Fn_mag = abs(Fn_c);
    F1_mag = abs(F1_c);
    th_hc  = -angle(Fn_c);
end

out.F_hc         = Fn_mag;
out.F_mc         = F1_mag;
out.theta_hc     = th_hc;
out.sigma_t_eq   = sig_t;
out.det_angle_hc = det_angle_hc;
out.psi_hc_deg   = psi_hc * 180.0 / pi;
out.U_hc_real    = 2.0 * (Fn_mag^2) * I0 * R_hc * cos(psi_hc)^2;
out.fais_mc_whc  = phis_cent;
out.V_hc_load_0  = -2.0 * Fn_mag * I0 * R_hc * cos(det_angle_hc) * exp(1i * (pi - det_angle_hc - th_hc));
out.xi_eff       = (2.0 * Fn_mag * I0 * R_hc * abs(cos(psi_hc))) / max(k_fp * V_mc, 1.0);
end
