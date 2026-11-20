import glob, re, scipy.io as sio, numpy as np

# Load all 290 files from the latest stopped scan
files = glob.glob('scan_data/Scan_Df50-250kHz_I300-500mA_20260921_105745/SOLEIL_*.mat')
print(f"Total files found: {len(files)}")

# Store classification results
old_codes = {}
new_codes = {}

for f in files:
    m = re.search(r'I0(\d+)mA.*detune(\d+)', f)
    if not m: continue
    i0_mA = int(m.group(1))
    df_Hz = int(m.group(2))
    key = (i0_mA, df_Hz)
    
    d = sio.loadmat(f)
    old_code = int(d['stability_code_point'][0,0]) if 'stability_code_point' in d else -1
    old_codes[key] = old_code
    
    qm = d['record_Q_mean'] # [3000, 416]
    pm = d['record_P_mean']
    qs = d['record_Q_std']
    
    sigma_t0 = float(d['sigma_t0'][0,0])
    sigma_e0 = float(d['sigma_e0'][0,0])
    Par_num = 10000
    I0 = i0_mA * 1e-3
    fre_shift_hc = df_Hz
    
    # Constants from HALF
    C = 353.98
    c = 299792458
    T0 = C / c
    tau_z = 13.7e-3
    damping_rate = 1 / tau_z
    
    nrec = qm.shape[0]
    win_steady = slice(int(nrec*0.7), nrec)
    win_len = nrec - int(nrec*0.7)
    
    # 1. All-bunch jitters
    bunch_jitters = np.std(qm[win_steady, :], axis=0) * sigma_t0 * 1e12 # ps
    mean_jit = np.mean(bunch_jitters)
    max_jit = np.max(bunch_jitters)
    
    # 2. Bunch-to-bunch pattern
    std_bb = np.std(np.mean(qm[-100:, :], axis=0)) * sigma_t0 * 1e12
    
    # 3. Bunch length oscillation
    bl_mean_t = np.mean(qs[win_steady, :], axis=1) * sigma_t0 * 1e12
    bl_osc_rel = np.std(bl_mean_t) / np.mean(bl_mean_t)
    
    # 4. Global all-bunch energy envelope
    P_detrend = pm[int(nrec*0.5):, :] - np.mean(pm[int(nrec*0.5):, :], axis=0)
    rms_P = np.sqrt(np.mean(P_detrend**2, axis=1)) # [win_len]
    
    # Noise floor
    noise_floor_ps = (sigma_t0 * 1e12) / np.sqrt(Par_num) # ~0.3 ps
    noise_floor_P  = 1.0 / np.sqrt(Par_num)              # ~0.01
    
    # Fit growth rate only if signal exceeds 2x noise floor
    pass_turns = np.arange(int(nrec*0.5), nrec) * 10
    if np.max(rms_P) > 2.5 * noise_floor_P:
        xx = np.log(rms_P + 1e-30)
        p_fit = np.polyfit(pass_turns, xx, 1)
        gr_si = p_fit[0] / T0
    else:
        gr_si = -10.0 # pure noise, strongly damped
        
    # 5. Coupled-bunch mode spatial FFT
    spatial_fft = np.fft.fft(qm[win_steady, :], axis=1)
    mode_power = np.mean(np.abs(spatial_fft)**2, axis=0)
    ac_modes = mode_power[1:]
    peak_mode = np.argmax(ac_modes) + 1
    crest_factor = ac_modes[peak_mode - 1] / np.mean(ac_modes)
    
    # 6. Physical classification
    R_hc = 1.86e6
    Q_hc = 31e3
    n_hc = 4
    U0 = 470e3
    V_mc = 1.7e6
    fre_hc = (416 / T0) * n_hc
    
    psi_hc = np.arctan(2 * Q_hc * fre_shift_hc / fre_hc)
    V_hc = 2 * I0 * R_hc * np.cos(psi_hc)
    k_fp = np.sqrt(1/n_hc**2 - 1/(n_hc**2-1) * (U0/V_mc)**2)
    xi_calc = V_hc / (k_fp * V_mc)
    
    # Beam is unstable if:
    # 1. Exponential growth exceeds damping
    # 2. OR sustained jitter significantly exceeds noise floor
    # 3. OR PTBL pattern is sustained
    # 4. OR coupled-bunch mode peak is prominent (crest factor > 5.0) AND max jitter > 0.6 ps
    has_lcbi_mode = (crest_factor > 5.0 and max_jit > 0.6) or (gr_si > 0.5 * damping_rate and max_jit > 0.5)
    
    is_unstable = (gr_si > 0.3 * damping_rate) or \
                  (mean_jit > max(1.2, 3.0 * noise_floor_ps)) or \
                  (max_jit > max(1.8, 5.0 * noise_floor_ps) and I0 >= 0.35) or \
                  (bl_osc_rel > 0.05) or \
                  (std_bb > 0.6) or \
                  has_lcbi_mode
                  
    code = 0
    if is_unstable:
        if (xi_calc >= 1.05 or (std_bb > 0.6 and np.mean(bl_mean_t) > 75)) and I0 >= 0.375:
            code = 1 # PTBL
        elif has_lcbi_mode and xi_calc < 0.98:
            code = 3 # LCBI
        elif bl_osc_rel > 0.05 or mean_jit > 1.0 or (xi_calc >= 0.85 and xi_calc <= 1.05):
            code = 2 # Coupled D-Q
        else:
            code = 4 # Other
            
    new_codes[key] = code

# Print summary comparison by current
print("\n" + "="*80)
print(f"{'I0 (mA)':>8} | {'OLD codes':<35} | {'NEW codes (Full-Bunch)':<35}")
print("="*80)

for i0_val in sorted(list(set(k[0] for k in old_codes))):
    keys = sorted([k for k in old_codes if k[0] == i0_val], key=lambda x: x[1])
    c_old = [old_codes[k] for k in keys]
    c_new = [new_codes[k] for k in keys]
    print(f"{i0_val:8d} | {str(c_old):<35} | {str(c_new):<35}")
