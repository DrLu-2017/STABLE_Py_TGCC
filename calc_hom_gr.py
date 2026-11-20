import numpy as np

# Machine parameters
E0 = 2.75e9
alpha_c = 1.07e-4
C = 353.98
c = 299792458
T0 = C / c
f0 = 1 / T0
w0 = 2 * np.pi * f0
h = 416
I0 = 0.3
tau_z = 13.7e-3  # damping time = 13.7 ms
gamma = E0 / 0.511e6
eta = alpha_c - 1/gamma**2

V_mc = 1.7e6
U0 = 470e3
Qs = np.sqrt(alpha_c * h / (2*np.pi) * np.sqrt((V_mc/E0)**2 - (U0/E0)**2))
ws = Qs * w0

print(f"f0 = {f0:.1f} Hz, Qs = {Qs:.4e}, fs = {Qs*f0:.1f} Hz")
print(f"Radiation Damping rate 1/tau_z = {1/tau_z:.1f} s^-1")

homs = [
    {'name': 'HOM 1', 'f': 1697.1e6, 'Q': 830, 'R': 6500},
    {'name': 'HOM 2', 'f': 1180.0e6, 'Q': 921, 'R': 7300},
    {'name': 'HOM 3', 'f': 959.1e6,  'Q': 209, 'R': 10100},
]

for hom in homs:
    fr = hom['f']
    wr = 2 * np.pi * fr
    Q = hom['Q']
    R = hom['R']
    
    max_gr = 0
    worst_mu = 0
    worst_val = 0
    
    # Mode mu from 0 to h-1
    for mu in range(h):
        # Revolution harmonics k = p*h + mu
        # Relevant frequencies are within a few bandwidths of fr
        BW = fr / Q
        f_min = fr - 5 * BW
        f_max = fr + 5 * BW
        k_min = int(np.floor(f_min / f0))
        k_max = int(np.ceil(f_max / f0))
        
        # Select k such that k mod h == mu
        gr_mu = 0
        for k in range(k_min, k_max + 1):
            if (k % h) == mu:
                w_plus = (k + Qs) * w0
                w_minus = (k - Qs) * w0
                
                # Resonator impedance Z(w)
                Z_p = R / (1 + 1j * Q * (w_plus/wr - wr/w_plus))
                Z_m = R / (1 + 1j * Q * (w_minus/wr - wr/w_minus))
                
                # Longitudinal dipole growth rate contribution:
                # Wang / Sacherer formula:
                # 1/tau = (eta * I0 / (2 * E0 * (Qs/w0))) * sum (w_p * Re(Z_p) - w_m * Re(Z_m))
                # Or Robinson:
                term = (w_plus * np.real(Z_p) - w_minus * np.real(Z_m))
                gr_mu += (eta * I0 / (2 * E0 * Qs * T0)) * term / (2*np.pi*f0) # normalized
                
        if gr_mu > max_gr:
            max_gr = gr_mu
            worst_mu = mu
            worst_val = gr_mu
            
    print(f"{hom['name']}: f={fr*1e-6:.1f} MHz, R={R} Ohm, Q={Q} -> worst mu={worst_mu}, growth_rate={worst_val:.1f} s^-1")
