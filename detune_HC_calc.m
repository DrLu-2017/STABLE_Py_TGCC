function [fre_detRQ] = detune_HC_calc(Ib, n_hc, C, h, U0, V_mc, R_hc, Q_hc, E0, alpha_c, sigma_e0, use_form_factor, xi_target)
% Compute the near-optimum HHC detuning frequency [Hz]
% for maximum bunch lengthening (flat-potential condition).
%
% Inputs:
%   Ib              - total beam current [A]
%   n_hc            - HHC harmonic order
%   C               - ring circumference [m]
%   h               - RF harmonic number
%   U0              - energy loss per turn [eV]
%   V_mc            - main cavity voltage [V]
%   R_hc            - HHC loaded shunt impedance [Ohm]
%   Q_hc            - HHC loaded quality factor
%   E0, alpha_c, sigma_e0 - (optional) beam energy [eV], mcf, natural energy spread
%   use_form_factor - (optional) bool, include iterative bunch form factor F_n(sigma_tau)
%   xi_target       - (optional) target flat-potential parameter (default 1.0)
%
% Output:
%   fre_detRQ - optimal HHC detuning frequency [Hz]  (positive = capacitive)

if nargin < 13 || isempty(xi_target)
    xi_target = 1.0;
end
if nargin < 12 || isempty(use_form_factor)
    use_form_factor = false;
end

if use_form_factor && nargin >= 11 && ~isempty(E0) && ~isempty(alpha_c) && ~isempty(sigma_e0)
    fp_out = calc_bunch_form_factor_eq('detune_for_xi', Ib, V_mc, U0, E0, C, h, alpha_c, sigma_e0, ...
        n_hc, R_hc, Q_hc, xi_target, true);
    fre_detRQ = round(fp_out.fre_shift_hc);
    disp([num2str(Ib*1e3), ' mA  ->  optimal HHC detuning (F_n=', num2str(fp_out.F_hc, '%.4f'), '): ', num2str(fre_detRQ), ' Hz']);
    return;
end

n2        = n_hc^2;
fre_hc    = 299792458 / C * h * n_hc;        % HHC resonant frequency [Hz]
k_fp      = sqrt(max(0, 1/n2 - 1/(n2-1) * (U0/V_mc)^2));
if Ib <= 0 || R_hc <= 0
    fre_detRQ = 0;
    return;
end
cos_arg   = (xi_target * k_fp * V_mc) / (-2*Ib*R_hc);
cos_arg   = max(-1, min(1, cos_arg));
phih_RQ   = acos(cos_arg);
fre_detRQ = round(tan(pi - phih_RQ) / (2*Q_hc) * fre_hc);
disp([num2str(Ib*1e3), ' mA  ->  optimal HHC detuning: ', num2str(fre_detRQ), ' Hz']);
end
