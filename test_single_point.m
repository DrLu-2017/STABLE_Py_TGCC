% Run single point tracking to diagnose instability and MC voltage
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
fre_shift_hc= 150e3;

% Check machine setup
cspeed = 299792458;
f0_est = cspeed / C;
Qs_est = sqrt(alpha_c / (2*pi) * h * sqrt((V_mc/E0)^2 - (U0/E0)^2));
T0_est = C / cspeed;
sigma_t0 = alpha_c * sigma_e0 / (2*pi*Qs_est) * T0_est;

fillrate = 1.0;
HALF = machine(C, I0, U0, E0, tau_s, tau_z, sigma_t0, sigma_e0, alpha_c, ...
               h, V_mc, n_hc, R_hc, Q_hc, fillrate, fre_shift_hc, ...
               Q_mc, R_mc, fre_shift_mc);

fprintf('\nHALF.Vrf_ideal  = %.3e V, angle = %.2f deg\n', abs(HALF.Vrf_ideal), angle(HALF.Vrf_ideal)*180/pi);
fprintf('HALF.V_mc_load_0= %.3e V, angle = %.2f deg\n', abs(HALF.V_mc_load_0), angle(HALF.V_mc_load_0)*180/pi);
fprintf('HALF.Vg_mc_init = %.3e V, angle = %.2f deg\n', abs(HALF.Vg_mc_init), angle(HALF.Vg_mc_init)*180/pi);

% Initial sum
Vc_init = HALF.Vg_mc_init + HALF.V_mc_load_0;
fprintf('Total Vc_init   = %.3e V (%.2f kV)\n', abs(Vc_init), abs(Vc_init)/1e3);

% Test PI settings
PI.KP        = 1;
PI.KI        = 1e-4;
PI.m         = 80;
PI.d         = 800;
PI.RL        = HALF.R_mc;
PI.Integral  = 0;
PI.DIg       = 0;
PI.Ig0       = HALF.Ig_mc_0;
PI.Ig        = HALF.Ig_mc_0;
PI.Vg0       = HALF.Vg_mc_0;
PI.DI        = 0;
PI.Ig_track  = [];
PI.Ig_FBid   = [];
PI.RoverQ    = HALF.R_mc / HALF.Q_mc;
PI.Vrf_mc_track = zeros(1, PI.m);
PI.enable    = true;
PI.Enabled   = true;

fprintf('PI.Ig0 = %.2f A, PI.Vg0 = %.3e V\n', abs(PI.Ig0), abs(PI.Vg0));
