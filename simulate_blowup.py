import numpy as np

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
fre_shift_hc = 150e3

f0 = cspeed / C
T0 = 1 / f0
Tb = T0 / h
f_rf = h * f0
w_rf = 2 * np.pi * f_rf
wr_mc = w_rf + 2 * np.pi * fre_shift_mc
wr_hc = n_hc * w_rf + 2 * np.pi * fre_shift_hc
wre_hc = n_hc * w_rf

Qs_est = np.sqrt(alpha_c / (2*np.pi) * h * np.sqrt((V_mc/E0)**2 - (U0/E0)**2))
sigma_t0 = alpha_c * sigma_e0 / (2*np.pi*Qs_est) * T0

det_angle_mc = np.arctan(Q_mc * (wr_mc/w_rf - w_rf/wr_mc))
det_angle_hc = np.pi - np.arctan(Q_hc * (wr_hc/wre_hc - wre_hc/wr_hc))

fais_mc_whc = np.pi - np.arcsin((U0 + 2*I0*R_hc*np.cos(det_angle_hc)**2) / V_mc)
Vrf_ideal = -V_mc * np.exp(-1j * (np.pi/2 - fais_mc_whc))
V_mc_load_0 = 2*I0*R_mc*np.cos(det_angle_mc)*np.exp(1j*det_angle_mc)
Vg_mc_init = Vrf_ideal - V_mc_load_0

Ig_mc_0 = Vg_mc_init / (2 * R_mc * np.cos(det_angle_mc)) * np.exp(-1j * det_angle_mc)
Vg_mc_0 = 2 * np.pi * (R_mc / Q_mc) * Ig_mc_0

V_hc_load_0 = -2*I0*R_hc*np.cos(det_angle_hc)*np.exp(1j*(np.pi - det_angle_hc))

drift_coef = alpha_c * T0 * (sigma_e0 / sigma_t0)
kick_coef  = 1.0 / (E0 * sigma_e0)
radampcoef = 2 * (T0 / tau_s)
quanexcoef = 2 * np.sqrt(T0 / tau_z)
ploss = kick_coef * U0

rot_decay_coef_hc = 1j*np.sqrt(1 - 1/(2*Q_hc)**2) - 1/(2*Q_hc)
TbAng_coef_hc = np.exp(rot_decay_coef_hc * (wr_hc * Tb))
exp_ang_coef_hc = -rot_decay_coef_hc * wr_hc * sigma_t0

rot_decay_coef_mc = 1j*np.sqrt(1 - 1/(2*Q_mc)**2) - 1/(2*Q_mc)
TbAng_coef_mc = np.exp(rot_decay_coef_mc * (wr_mc * Tb))
exp_ang_coef_mc = -rot_decay_coef_mc * wr_mc * sigma_t0

VbImagFactor_hc = 1/np.sqrt(4*Q_hc**2 - 1)
VbImagFactor_mc = 1/np.sqrt(4*Q_mc**2 - 1)

Par_num = 200 # fast test
Bun_num = 416
qc = T0 * I0 / h / Par_num
Vb_hc = qc * wr_hc * R_hc / Q_hc
Vb_mc = qc * wr_mc * R_mc / Q_mc

# Initial distribution
np.random.seed(42)
Q = np.random.randn(Par_num, Bun_num)
P = np.random.randn(Par_num, Bun_num)

# PI setup
PI_m = 80
PI_d = 800
PI_KP = 1.0
PI_KI = 1e-4
PI_RL = R_mc
PI_RoverQ = R_mc / Q_mc
piIndex = 0
PI_Integral = 0.0
PI_Vg0 = Vg_mc_0
PI_Ig = Ig_mc_0
PI_Ig_FBid = []
PI_Vrf_mc_track = np.zeros(PI_m, dtype=complex)

Vg_mc_track_0 = Vg_mc_init
V_hc_load = V_hc_load_0
V_mc_load = V_mc_load_0

print("Starting 1000 turn test...")
for turn in range(1, 1001):
    Q += P * drift_coef
    
    # HHC
    exp_angle = np.exp(exp_ang_coef_hc * Q)
    exp_sum = np.sum(exp_angle, axis=0) # [416]
    V_load_cpu = exp_sum * Vb_hc
    
    # Cavity rotation through turn
    V_hc_load_turn = np.zeros(Bun_num, dtype=complex)
    for b in range(Bun_num):
        V_hc_load = V_hc_load * TbAng_coef_hc + V_load_cpu[b]
        V_hc_load_turn[b] = V_hc_load
    V_hc_kick = (V_hc_load_turn * kick_coef) / exp_angle
    
    # MC
    exp_angle_mc = np.exp(exp_ang_coef_mc * Q)
    exp_sum_mc = np.sum(exp_angle_mc, axis=0)
    V_load_mc_cpu = exp_sum_mc * Vb_mc
    
    Vc_mc = np.zeros(Bun_num, dtype=complex)
    for b in range(Bun_num):
        piIndex += 1
        if piIndex > PI_d and (piIndex - PI_d) % PI_m == 0:
            if len(PI_Ig_FBid) > 0:
                PI_Vg0 = 2 * np.pi * PI_RoverQ * PI_Ig_FBid.pop(0)
        
        Vg_mc_track_0 = Vg_mc_track_0 * TbAng_coef_mc + PI_Vg0
        V_mc_load = V_mc_load * TbAng_coef_mc + V_load_mc_cpu[b]
        Vc_mc[b] = Vg_mc_track_0 + V_mc_load
        
        piTrackIndex = (piIndex - 1) % PI_m
        PI_Vrf_mc_track[piTrackIndex] = Vc_mc[b]
        if piTrackIndex == PI_m - 1:
            Vrf_mean = np.mean(PI_Vrf_mc_track)
            DI = (Vrf_ideal - Vrf_mean) / PI_RL
            PI_Integral += PI_KI * DI
            DIg = PI_KP * DI + PI_Integral
            PI_Ig = Ig_mc_0 + DIg
            PI_Ig_FBid.append(PI_Ig)
            
    V_mc_kick = (Vc_mc * kick_coef) / exp_angle_mc
    
    # Radiation damping & kick
    rad_kick = -radampcoef * P + quanexcoef * np.random.randn(Par_num, Bun_num)
    P = P + rad_kick - ploss - np.real(V_hc_kick) + np.imag(V_hc_kick)*VbImagFactor_hc \
        - np.real(V_mc_kick) + np.imag(V_mc_kick)*VbImagFactor_mc
        
    mean_tau_ps = np.mean(Q[:, 0]) * sigma_t0 * 1e12
    std_tau_ps  = np.std(Q[:, 0]) * sigma_t0 * 1e12
    vc_amp_kV   = np.abs(Vc_mc[0]) / 1e3
    
    if turn in [1, 10, 50, 100, 200, 400, 600, 700, 750, 800, 820, 840, 850]:
        print(f"Turn {turn:4d}: <tau> = {mean_tau_ps:8.2f} ps, std_tau = {std_tau_ps:6.2f} ps, |Vc_mc| = {vc_amp_kV:10.2f} kV")
        if np.abs(mean_tau_ps) > 5000:
            print(f"--> BEAM BLOWN UP AT TURN {turn}!")
            break
