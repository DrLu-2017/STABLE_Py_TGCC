% =========================================================================
% machine_Soleil2.m
% Machine parameters for Soleil 2 (4th generation storage ring, France)
% Compatible with STABLE tracking code (MPMBTrackBunchLengthingwoHighQ.m)
%
% Usage: replace the beam parameter block at the top of
%        MPMBTrackBunchLengthingwoHighQ.m with this script,
%        OR call it as a standalone setup before the tracking loop.
%
% Reference: Soleil 2 CDR / lattice design reports (2021-2024)
% Last updated: 2026-06
% =========================================================================

%% -------------------------------------------------------------------------
%% 1. Fundamental beam / lattice parameters
%% -------------------------------------------------------------------------
cspeed   = 299792458;       % speed of light [m/s]

C        = 354.097;         % [m]      ring circumference
E0       = 2.75e9;          % [eV]     beam energy  (2.75 GeV)
I0       = 500e-3;          % [A]      total stored beam current (design)
U0       = 456e3;           % [eV]     energy loss per turn (synchrotron radiation)
h        = 416;             % [ ]      RF harmonic number  (f_rf = h * f0)
alpha_c  = 1.27e-4;         % [ ]      momentum compaction factor

%% Radiation damping & equilibrium emittances
tau_s    = 7.3e-3;          % [s]      longitudinal radiation damping time
tau_z    = 7.3e-3;          % [s]      energy  radiation damping time
%   (for a 4GSR with strong damping wigglers tau can be ~3-5 ms;
%    7.3 ms corresponds to the bare lattice without wigglers)

sigma_e0 = 1.06e-3;         % [ ]      natural rms energy spread
%   sigma_t0 (natural bunch length without HHC):
%   sigma_t0 = alpha_c * sigma_e0 / omega_s * T0
%   Pre-compute here for reference; the value used in tracking is below.
f0_est   = cspeed / C;
f_rf_est = h * f0_est;
V_mc_check = 1.8e6;         % [V]  (used only for the sigma_t0 estimate below)
Qs_est   = sqrt( alpha_c / (2*pi) * h * sqrt((V_mc_check/E0)^2 - (U0/E0)^2) );
T0_est   = C / cspeed;
sigma_t0 = alpha_c * sigma_e0 / (2*pi*Qs_est) * T0_est;  % [s]  auto-computed
fprintf('Natural bunch length (no HHC): sigma_t0 = %.2f ps\n', sigma_t0*1e12);

%% -------------------------------------------------------------------------
%% 2. Main RF cavity (MC) — superconducting (SRF)
%% -------------------------------------------------------------------------
V_mc          = 1.8e6;      % [V]    total main cavity voltage
Q_mc          = 2.0e5;      % [ ]    loaded Q  (SRF, typical beta_c ~ 2)
R_over_Q_mc   = 45;         % [Ohm]  R/Q of main cavity (geometry factor)
R_mc          = Q_mc * R_over_Q_mc;  % [Ohm]  loaded shunt impedance

% Detuning of MC (inductive detuning to compensate beam loading):
%   A first estimate is: Df_mc = -Ib*R_mc / (2*Q_mc*V_mc) * f_rf  [Hz]
%   Fine-tune this value after running the first simulation.
fre_shift_mc  = -8.0e3;     % [Hz]   MC detuning (negative = inductive)

%% -------------------------------------------------------------------------
%% 3. Harmonic Cavity (HHC) — passive 3rd harmonic
%% -------------------------------------------------------------------------
n_hc     = 3;               % [ ]    harmonic order (3rd harmonic)
Q_hc     = 2.0e8;           % [ ]    loaded Q  (high-Q passive cavity)
R_over_Q_hc = 45;           % [Ohm]  R/Q of HHC (geometry factor)
R_hc     = Q_hc * R_over_Q_hc;  % [Ohm]  loaded shunt impedance

% Optimal HHC detuning — computed from near-flat-potential condition:
fre_shift_hc = detune_HC_calc(I0, n_hc, C, h, U0, V_mc, R_hc, Q_hc);
%   If you want to scan detuning instead, comment the line above and use:
%   fre_shift_scan = [60:20:200] * 1e3;   % [Hz]
%   fre_shift_hc   = fre_shift_scan(Ii);  % inside the scan loop

%% -------------------------------------------------------------------------
%% 4. Fill pattern
%% -------------------------------------------------------------------------
% --- Option A: Uniform filling (simplest, good starting point) ---
pattern  = ones(1, h);

% --- Option B: ~80% uniform filling (8 empty buckets per 40-bucket block) ---
% pattern  = ones(1, h);
% pat = 33:40;
% for i = 2 : h/40
%     pat = [pat, 40*i-7 : 40*i];
% end
% pattern(pat) = 0;

% --- Option C: Single bunch (for growth rate tests) ---
% pattern      = zeros(1, h);
% pattern(1)   = 1;

fillrate = sum(pattern) / h;
fprintf('Fill rate: %.1f %%  (%d bunches out of %d buckets)\n', ...
        fillrate*100, sum(pattern), h);

%% -------------------------------------------------------------------------
%% 5. Call machine() to build the HALF structure
%% -------------------------------------------------------------------------
HALF = machine(C, I0, U0, E0, tau_s, tau_z, ...
               sigma_t0, sigma_e0, alpha_c, ...
               h, V_mc, ...
               n_hc, R_hc, Q_hc, ...
               fillrate, fre_shift_hc, ...
               Q_mc, R_mc, fre_shift_mc);

%% -------------------------------------------------------------------------
%% 6. Short-range wake switch
%% -------------------------------------------------------------------------
HALF.ShortRange_on = 1;   % 1 = include intra-bunch short-range wake
                           % 0 = disable (faster, use for initial tests)

%% -------------------------------------------------------------------------
%% 7. Sanity-check printout
%% -------------------------------------------------------------------------
fprintf('\n===== Soleil 2 — machine() output summary =====\n');
fprintf('Ring circumference     C   = %.3f m\n',   HALF.C);
fprintf('Revolution period      T0  = %.4f us\n',  HALF.T0*1e6);
fprintf('Revolution freq        f0  = %.4f kHz\n', 1/HALF.T0*1e-3);
fprintf('RF frequency           frf = %.4f MHz\n', HALF.f_rf*1e-6);
fprintf('HHC frequency          fhc = %.4f MHz\n', HALF.fre_hc*1e-6);
fprintf('HHC detuning           Df  = %.2f kHz\n', fre_shift_hc*1e-3);
fprintf('Bunch charge           qc  = %.4f nC\n',  HALF.qc*1e9);
fprintf('Natural synch. phase   phi = %.3f deg\n', HALF.fais_nat/pi*180);
fprintf('MC synch. phase (w/HHC)    = %.3f deg\n', HALF.fais_mc_whc/pi*180);
fprintf('MC detuning angle          = %.4f deg\n', HALF.det_angle_mc/pi*180);
fprintf('HHC detuning angle         = %.4f deg\n', HALF.det_angle_hc/pi*180);
fprintf('Natural bunch length sigma_t0 = %.2f ps\n', sigma_t0*1e12);
fprintf('================================================\n\n');
