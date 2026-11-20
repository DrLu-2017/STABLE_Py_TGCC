import numpy as np

# Verify that open-loop generator drive maintains cavity voltage and stable beam
E0 = 2.75e9
U0 = 470e3
sigma_e0 = 0.92e-3
C = 353.98
c = 299792458
T0 = C / c
h = 416
Tb = T0 / h
f_rf = h / T0
w_rf = 2 * np.pi * f_rf
alpha_c = 1.07e-4

V_mc = 1.7e6
Q_mc = 35.7e3
RoQ_mc = 45.0
R_mc = Q_mc * RoQ_mc
fre_shift_mc = -35.0e3
wre_mc = w_rf
wr_mc = wre_mc + 2 * np.pi * fre_shift_mc

I0 = 0.3 # 300 mA
fais_nat = np.pi - np.arcsin(U0 / V_mc)

# Detuning angle
det_angle_mc = np.arctan(Q_mc * (wr_mc/wre_mc - wre_mc/wr_mc))
# Beam loading voltage
V_mc_load_0 = 2 * I0 * R_mc * np.cos(det_angle_mc) * np.exp(1j * det_angle_mc)

Vrf_ideal = -V_mc * np.exp(-1j * (np.pi/2 - fais_nat))
Vg_mc_init = Vrf_ideal - V_mc_load_0
Ig_mc_0 = Vg_mc_init / (2 * R_mc * np.cos(det_angle_mc)) * np.exp(-1j * det_angle_mc)
Vg_mc_0 = 2 * np.pi * (R_mc / Q_mc) * Ig_mc_0

rot_decay_coef_mc = 1j * np.sqrt(1 - 1/(2*Q_mc)**2) - 1/(2*Q_mc)
TbAng_coef_mc = np.exp(rot_decay_coef_mc * wr_mc * Tb)

print("--- Test 1: With Vg0 = 0 (the BUG when use_PI=false) ---")
Vg = Vg_mc_init
for turn in range(2000):
    for bucket in range(h):
        Vg = Vg * TbAng_coef_mc + 0  # BUG!
    if turn in [0, 10, 50, 100, 500, 1000, 1800]:
        print(f"Turn {turn:4d}: |Vg| = {abs(Vg):12.2f} V (angle = {np.angle(Vg)*180/np.pi:6.1f} deg)")

print("\n--- Test 2: With fixed feedforward Vg0 = Vg_mc_0 (the FIX) ---")
Vg = Vg_mc_init
for turn in range(2000):
    for bucket in range(h):
        Vg = Vg * TbAng_coef_mc + Vg_mc_0
    if turn in [0, 10, 50, 100, 500, 1000, 1800]:
        print(f"Turn {turn:4d}: |Vg| = {abs(Vg):12.2f} V (angle = {np.angle(Vg)*180/np.pi:6.1f} deg)")
