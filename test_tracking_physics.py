import numpy as np

# Reproduce STABLE tracking in Python to diagnose why single points are unstable
cspeed = 299792458
C = 354.097
E0 = 2.75e9
I0 = 0.5
U0 = 470e3
h = 416
alpha_c = 1.07e-4
tau_s = 11.6e-3
tau_z = 13.7e-3
sigma_e0 = 1.06e-3

V_mc = 1.7e6
Q_mc = 35.7e3
RoQ_mc = 140
R_mc = 5e6
fre_shift_mc = -35e3

n_hc = 4
Q_hc = 31e3
RoQ_hc = 30
R_hc = 1.86e6

f0 = cspeed / C
T0 = 1 / f0
Tb = T0 / h
f_rf = h * f0
w_rf = 2 * np.pi * f_rf
wr_mc = w_rf + 2 * np.pi * fre_shift_mc
wre_hc = n_hc * w_rf

Qs_est = np.sqrt(alpha_c / (2*np.pi) * h * np.sqrt((V_mc/E0)**2 - (U0/E0)**2))
sigma_t0 = alpha_c * sigma_e0 / (2*np.pi*Qs_est) * T0

det_angle_mc = np.arctan(Q_mc * (wr_mc/w_rf - w_rf/wr_mc))

# Test two detunings: 100 kHz (flat-pot) and 150 kHz
for fre_shift_hc in [100e3, 150e3]:
    wr_hc = n_hc * w_rf + 2 * np.pi * fre_shift_hc
    det_angle_hc = np.pi - np.arctan(Q_hc * (wr_hc/wre_hc - wre_hc/wr_hc))
    
    fais_mc_whc = np.pi - np.arcsin((U0 + 2*I0*R_hc*np.cos(det_angle_hc)**2) / V_mc)
    Vrf_ideal = -V_mc * np.exp(-1j * (np.pi/2 - fais_mc_whc))
    V_mc_load_0 = 2*I0*R_mc*np.cos(det_angle_mc)*np.exp(1j*det_angle_mc)
    Vg_mc_init = Vrf_ideal - V_mc_load_0
    
    Ig_mc_0 = Vg_mc_init / (2 * R_mc * np.cos(det_angle_mc)) * np.exp(-1j * det_angle_mc)
    Vg_mc_0 = 2 * np.pi * (R_mc / Q_mc) * Ig_mc_0
    
    V_hc_load_0 = -2*I0*R_hc*np.cos(det_angle_hc)*np.exp(1j*(np.pi - det_angle_hc))
    
    print(f"\n--- Testing df_hc = {fre_shift_hc/1e3:.1f} kHz ---")
    print(f"|Vrf_ideal| = {np.abs(Vrf_ideal):.3e}, |V_mc_load_0| = {np.abs(V_mc_load_0):.3e}, |Vg_mc_init| = {np.abs(Vg_mc_init):.3e}")
    print(f"det_angle_mc = {det_angle_mc*180/np.pi:.2f} deg, det_angle_hc = {det_angle_hc*180/np.pi:.2f} deg")
