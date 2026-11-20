% Test HOM sideband alignment formula
cspeed_val = 299792458;
C_val      = 353.98;
h_val      = 416;
E0_val     = 2.75e9;
alpha_val  = 1.07e-4;
V_mc_val   = 1.7e6;
U0_val     = 470e3;

f0_val = cspeed_val / C_val;
sin_phi_s = min(1, max(-1, U0_val / V_mc_val));
cos_phi_s = sqrt(max(0, 1 - sin_phi_s^2));
nu_s = sqrt((alpha_val * h_val * V_mc_val * cos_phi_s) / (2 * pi * E0_val));
fs_val = nu_s * f0_val;

fprintf('f0 = %.2f Hz, nu_s = %.5f, fs = %.2f Hz\n', f0_val, nu_s, fs_val);

f_noms = [1697.1e6; 1180.0e6; 959.1e6];
for idx = 1:length(f_noms)
    fn = f_noms(idx);
    k = round((fn - fs_val) / f0_val);
    M = floor(k / h_val);
    mu = mod(k, h_val);
    f_align = (k + nu_s) * f0_val;
    fprintf('Nominal: %.1f MHz -> k=%d (M=%d, mu=%d), f_align = %.3f Hz (shift = %+.1f kHz)\n', ...
        fn/1e6, k, M, mu, f_align, (f_align - fn)/1e3);
end
