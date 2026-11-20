function [HALF]=machine(C,I0,U0,E0,tau_s,tau_z,sigma_t0,sigma_e0,alpha_c,h,V_mc,n_hc,R_hc,Q_hc,fillrate,fre_shift_hc,Q_mc,R_mc,fre_shift_mc,use_form_factor)
% Build the HALF machine structure with all derived physical quantities.
% All disp() messages in English to avoid encoding issues.
if nargin < 20 || isempty(use_form_factor)
    use_form_factor = false;
end

cspeed = 299792458;
HALF.C       = C;
HALF.R       = C/(2*pi);
HALF.I0      = I0;
HALF.h       = h;
HALF.U0      = U0;
HALF.E0      = E0;
HALF.sigma_t0= sigma_t0;
HALF.sigma_e0= sigma_e0;
HALF.alpha_c = alpha_c;

% Enforce longitudinal Fokker-Planck fluctuation-dissipation consistency (<P^2>_eq = tau_s / tau_z == 1).
% If distinct values are supplied, use the worst-case (longer, weaker-damping) longitudinal damping time.
if abs(tau_s - tau_z) > 1e-9
    tau_long = max(tau_s, tau_z);
    disp(['[Machine] Notice: tau_s (', num2str(tau_s*1e3), ' ms) != tau_z (', num2str(tau_z*1e3), ...
          ' ms); enforcing Fokker-Planck balance with worst-case tau_z = ', num2str(tau_long*1e3), ' ms']);
    tau_s = tau_long;
    tau_z = tau_long;
end
HALF.tau_s   = tau_s;
HALF.tau_z   = tau_z;
HALF.V_mc    = V_mc;
HALF.n_hc    = n_hc;
HALF.R_hc    = R_hc;
HALF.Q_hc    = Q_hc;
HALF.fillrate  = fillrate;
HALF.fre_shift_hc = fre_shift_hc;
HALF.R_mc    = R_mc;
HALF.Q_mc    = Q_mc;
HALF.fre_shift_mc = fre_shift_mc;
HALF.use_form_factor = logical(use_form_factor);

HALF.gamma = E0/0.51099906/1e6;       % Lorentz gamma  (electron rest mass = 0.511 MeV)
HALF.beta  = sqrt(1-1/(HALF.gamma^2));

HALF.T0   = HALF.C/cspeed/HALF.beta;  % revolution period [s]
HALF.Tb   = HALF.T0/HALF.h;           % bucket spacing [s]
HALF.f_rf = HALF.h/HALF.T0;           % RF frequency [Hz]
HALF.w_rf = HALF.f_rf*2*pi;           % RF angular frequency [rad/s]

HALF.fre_hc = HALF.f_rf * HALF.n_hc;  % HHC resonant frequency [Hz]
HALF.wre_hc = HALF.fre_hc * 2*pi;     % HHC resonant angular frequency [rad/s]
HALF.qc     = HALF.T0*HALF.I0/HALF.h/HALF.fillrate;  % bunch charge [C]
HALF.wr_hc  = HALF.w_rf*HALF.n_hc + fre_shift_hc*pi*2; % HHC actual angular frequency
HALF.angle_hc = HALF.wr_hc * HALF.Tb; % HHC phase advance per bucket [rad]

HALF.fre_mc   = HALF.f_rf;            % MC resonant frequency [Hz]
HALF.wre_mc   = HALF.fre_mc * 2*pi;   % MC resonant angular frequency [rad/s]
HALF.wr_mc    = HALF.wre_mc + fre_shift_mc*pi*2; % MC actual angular frequency
HALF.angle_mc = HALF.wr_mc * HALF.Tb; % MC phase advance per bucket [rad]

HALF.lambda_rads = HALF.T0 / HALF.tau_s;  % longitudinal radiation damping rate per turn
HALF.lambda_radz = HALF.T0 / HALF.tau_z;  % energy radiation damping rate per turn

HALF.fais_nat = pi - asin(HALF.U0/HALF.V_mc); % natural synchronous phase [rad]

%% --- HHC beam loading (point-bunch or iterative self-consistent form factor F_n) ---
HALF.det_angle_hc = pi - atan(HALF.Q_hc*(HALF.wr_hc/HALF.wre_hc - HALF.wre_hc/HALF.wr_hc));
disp(['[HHC] Detuning angle                  : ', num2str(HALF.det_angle_hc/pi*180), ' deg']);

if HALF.use_form_factor
    fp_eq = calc_bunch_form_factor_eq('eval', HALF.I0, HALF.V_mc, HALF.U0, HALF.E0, ...
        HALF.C, HALF.h, HALF.alpha_c, HALF.sigma_e0, HALF.n_hc, HALF.R_hc, HALF.Q_hc, HALF.fre_shift_hc, true);
    HALF.F_hc_eq    = fp_eq.F_hc;
    HALF.F_mc_eq    = fp_eq.F_mc;
    HALF.sigma_t_eq = fp_eq.sigma_t_eq;
    disp(['[HHC] Iterative bunch form factor F_n : ', num2str(HALF.F_hc_eq, '%.4f'), ...
          '  (sigma_t,eq = ', num2str(HALF.sigma_t_eq*1e12, '%.2f'), ' ps)']);
else
    HALF.F_hc_eq    = 1.0;
    HALF.F_mc_eq    = 1.0;
    HALF.sigma_t_eq = HALF.sigma_t0;
end
% Note: At turn 0, particles are initialized with natural bunch length sigma_t0 (prior to adiabatic
% lengthening), AND mbtrack2 sets MC.theta = pi - asin((U0 + 2*I0*R_hc*cos^2(psi_hc)) / V_mc).
% Using point-bunch initial cavity phasor & MC reference phase avoids turn-0 dipole shock while
% keeping F_hc_eq / sigma_t_eq for self-consistent flat-potential xi_eff & detuning calculations.
HALF.fais_mc_whc = pi - asin((HALF.U0 + 2*HALF.I0*HALF.R_hc*cos(HALF.det_angle_hc)^2) / HALF.V_mc);
HALF.V_hc_load_0 = -2*HALF.I0*HALF.R_hc*cos(HALF.det_angle_hc)*exp(1i*(pi - HALF.det_angle_hc));
disp(['[HHC] MC sync. phase (with HHC load)  : ', num2str(HALF.fais_mc_whc), ' rad  = ', num2str(HALF.fais_mc_whc/pi*180), ' deg']);
disp(['[HHC] Initial beam-loading voltage     : ', num2str(HALF.V_hc_load_0), ' V']);

%% --- MC beam loading ---
HALF.det_angle_mc = atan(HALF.Q_mc*(HALF.wr_mc/HALF.wre_mc - HALF.wre_mc/HALF.wr_mc));
disp(['[MC]  Detuning angle                   : ', num2str(HALF.det_angle_mc/pi*180), ' deg']);

HALF.V_mc_load_0 = 2*HALF.I0*HALF.R_mc*cos(HALF.det_angle_mc)*exp(1i*HALF.det_angle_mc);
disp(['[MC]  Initial beam-loading voltage     : ', num2str(HALF.V_mc_load_0), ' V']);

%% --- Generator voltage initialisation ---
HALF.Vrf_ideal      = -HALF.V_mc * exp(-1i*(pi/2 - HALF.fais_mc_whc));
HALF.Vg_mc_init     = HALF.Vrf_ideal - HALF.V_mc_load_0;
HALF.Vg_mc_track_0  = HALF.Vg_mc_init; % initial generator voltage phasor

HALF.Ig_mc_0 = HALF.Vg_mc_track_0 / (2*HALF.R_mc*cos(HALF.det_angle_mc)) * exp(-1i*HALF.det_angle_mc);
HALF.Vg_mc_0 = 2*pi*HALF.R_mc/HALF.Q_mc * HALF.Ig_mc_0;

%% --- Normalised tracking coefficients ---
HALF.drift_coef = HALF.alpha_c * HALF.T0 * (HALF.sigma_e0 / HALF.sigma_t0);
HALF.kick_coef  = 1 / (HALF.E0 * HALF.sigma_e0);

HALF.rfcoef1 = HALF.kick_coef * HALF.V_mc; % scaled MC voltage kick coefficient
HALF.rfcoef2 = HALF.w_rf * HALF.sigma_t0;  % normalised RF phase argument

HALF.radampcoef = 2 * HALF.lambda_rads;    % radiation damping coefficient
HALF.quanexcoef = 2 * sqrt(HALF.lambda_radz); % quantum excitation coefficient

HALF.ploss = HALF.kick_coef * HALF.U0;    % normalised energy loss per turn

%% --- High-Q cavity rotation / decay coefficients ---
HALF.VbImagFactor_hc = 1/sqrt(4*Q_hc^2 - 1);
HALF.rot_coef_hc     = sqrt(1 - 1/(2*Q_hc)^2);

HALF.VbImagFactor_mc = 1/sqrt(4*Q_mc^2 - 1);
HALF.rot_coef_mc     = sqrt(1 - 1/(2*Q_mc)^2);

%% --- Loading angle ---
HALF.angle_V_mc = Vb_angle_calc(real(HALF.Vrf_ideal), imag(HALF.Vrf_ideal));
HALF.angle_V_ge = Vb_angle_calc(real(HALF.Vg_mc_init), imag(HALF.Vg_mc_init));

HALF.angle_load = HALF.fais_mc_whc - (pi/2 - HALF.angle_V_ge - HALF.det_angle_mc);
disp(['[MC]  Initial loading angle            : ', num2str(HALF.angle_load/pi*180), ' deg']);
end
