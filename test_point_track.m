% Test single point tracking under various conditions
clear; clc;

% Load Soleil2 parameters
C        = 354.097;
E0       = 2.75e9;
I0       = 0.5;
U0       = 470e3;
h        = 416;
alpha_c  = 1.07e-4;
tau_s    = 11.6e-3;
tau_z    = 13.7e-3;
sigma_e0 = 1.06e-3;

V_mc        = 1.7e6;
Q_mc        = 35.7e3;
RoQ_mc      = 140;
R_mc        = 5e6;
fre_shift_mc= -35e3;

n_hc        = 4;
Q_hc        = 31e3;
RoQ_hc      = 30;
R_hc        = 1.86e6;

cspeed = 299792458;
Qs_est = sqrt(alpha_c / (2*pi) * h * sqrt((V_mc/E0)^2 - (U0/E0)^2));
T0_est = C / cspeed;
sigma_t0 = alpha_c * sigma_e0 / (2*pi*Qs_est) * T0_est;
fillrate = 1.0;
pattern = ones(1, h);

% Test Case 1: df_hc = 100 kHz (flat-pot xi~1), DFB ON (gain=0.2, every=1, delay=704)
% Test Case 2: df_hc = 150 kHz (user screenshot), DFB ON
% Test Case 3: df_hc = 150 kHz, DFB OFF

fprintf('=== Diagnostic Test: Checking initial voltages and Vc turn 1 ===\n');
for df = [100e3, 150e3]
    HALF = machine(C, I0, U0, E0, tau_s, tau_z, sigma_t0, sigma_e0, alpha_c, ...
                   h, V_mc, n_hc, R_hc, Q_hc, fillrate, df, ...
                   Q_mc, R_mc, fre_shift_mc);
    Vc_init_kV = abs(HALF.Vg_mc_init + HALF.V_mc_load_0) / 1e3;
    fprintf('df_hc = %6.1f kHz: Vc_init = %8.1f kV  (Target = 1700 kV)\n', df/1e3, Vc_init_kV);
end
