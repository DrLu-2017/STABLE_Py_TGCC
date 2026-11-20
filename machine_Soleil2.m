% =========================================================================
% machine_Soleil2.m (Worst-case RF stability configuration: IDs open)
% =========================================================================

%% 1. Fundamental beam / lattice parameters (Worst-case for RF: IDs open)
cspeed   = 299792458;
C        = 354.097;         % SOLEIL II nominal circumference [m] (f_rf = 352.20 MHz)
E0       = 2.75e9;
I0       = 500e-3;
U0       = 470e3;           % Bare lattice / IDs open energy loss per turn [eV]
h        = 416;
alpha_c  = 1.07e-4;

% Unified longitudinal radiation damping time (IDs open = 13.7 ms, weakest radiation damping)
tau_s    = 13.7e-3;
tau_z    = 13.7e-3;
sigma_e0 = 0.92e-3;         % IDs open natural RMS energy spread (weakest Landau damping)

f0_est   = cspeed / C;
f_rf_est = h * f0_est;
V_mc_check = 1.7e6;
Qs_est   = sqrt( alpha_c / (2*pi) * h * sqrt((V_mc_check/E0)^2 - (U0/E0)^2) );
T0_est   = C / cspeed;
sigma_t0 = alpha_c * sigma_e0 / (2*pi*Qs_est) * T0_est;

%% 2. Main RF cavity (MC)
V_mc          = 1.7e6;      % 1.7 MV
Ncav_mc       = 4;          % 4 main cavities
Q0_mc         = 35.7e3;     % Unloaded Q0
R_over_Q_mc   = 140;        % Ohm per cavity
beta_mc       = 4.95;       % Power coupler beta
Q_mc          = Q0_mc / (1 + beta_mc); % Loaded Q_L = 6000
R_mc          = Ncav_mc * Q_mc * R_over_Q_mc; % 3.36 MOhm total loaded impedance
fre_shift_mc  = -35.0e3;    % -35 kHz

%% 3. Harmonic Cavity (HHC)
n_hc        = 4;            % 4th harmonic (1.4088 GHz)
Ncav_hc     = 2;            % 2 passive cavities
Q0_hc       = 31e3;         % Unloaded Q0 (warm copper)
R_over_Q_hc = 30;           % Ohm per cavity
beta_hc     = 0;            % Passive cavity (no power coupler, beta=0)
Q_hc        = Q0_hc / (1 + beta_hc); % Loaded Q_L = 31000
R_hc        = Ncav_hc * Q_hc * R_over_Q_hc; % 1.86 MOhm total loaded impedance
fre_shift_hc = 95e3;        % ~95 kHz (near xi_eff ~ 1.0 at 500 mA with F_4)
use_form_factor = true;     % Iterative self-consistent bunch form factor F_n(sigma_tau)

%% 4. Fill pattern
pattern  = ones(1, h);
fillrate = sum(pattern) / h;

%% 5. Call machine()
SOLEILII = machine(C, I0, U0, E0, tau_s, tau_z, ...
               sigma_t0, sigma_e0, alpha_c, ...
               h, V_mc, ...
               n_hc, R_hc, Q_hc, ...
               fillrate, fre_shift_hc, ...
               Q_mc, R_mc, fre_shift_mc, use_form_factor);

%% 6. Short-range wake switch
HALF = SOLEILII;
HALF.ShortRange_on = 1;

%% 7. Sanity-check printout
fprintf('\n===== Soleil II machine() output summary =====\n');
fprintf('RF frequency           frf = %.4f MHz\n', (h * cspeed / C) * 1e-6);
fprintf('MC Voltage             Vc  = %.2f MV\n', V_mc*1e-6);
fprintf('MC Detuning            Df  = %.2f kHz\n', fre_shift_mc*1e-3);
fprintf('HHC Detuning           Df  = %.2f kHz (F_4 = %.4f)\n', fre_shift_hc*1e-3, SOLEILII.F_hc_eq);
fprintf('================================================\n\n');
