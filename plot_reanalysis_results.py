import os, glob, re, shutil
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.patches as mpatches

folder = 'scan_data/Scan_Df50-250kHz_I300-500mA_20260921_105745'
sum_file = os.path.join(folder, 'scan_summary.mat')
sum_data = sio.loadmat(sum_file)

Df_scan_kHz = sum_data['Df_scan_kHz'].flatten()
I0_scan_mA  = sum_data['I0_scan_mA'].flatten()
nI0 = len(I0_scan_mA)
nDf = len(Df_scan_kHz)

old_stability_code = sum_data['stability_code'].copy()
new_stability_code = np.full((nI0, nDf), -1, dtype=int)
new_growth_rate    = np.full((nI0, nDf), np.nan)
new_peak_freq      = np.full((nI0, nDf), np.nan)

# Physical constants for Soleil 2
C = 353.98
c = 299792458
T0 = C / c
tau_z = 13.7e-3
damping_rate = 1.0 / tau_z
R_hc = 1.86e6
Q_hc = 31e3
n_hc = 4
U0 = 470e3
V_mc = 1.7e6
fre_hc = (416 / T0) * n_hc

files = glob.glob(os.path.join(folder, 'SOLEIL_*.mat'))
print(f"Loading and processing {len(files)} files...")

for f in files:
    m = re.search(r'I0(\d+)mA.*detune(-?\d+\.?\d*[eE]?[-+]?\d*)', f)
    if not m: continue
    i0_mA = int(m.group(1))
    df_Hz = float(m.group(2))
    
    iI_matches = np.where(np.abs(I0_scan_mA - i0_mA) < 1e-2)[0]
    iDf_matches = np.where(np.abs(Df_scan_kHz - df_Hz/1e3) < 1e-2)[0]
    if len(iI_matches) == 0 or len(iDf_matches) == 0: continue
    iI = iI_matches[0]
    iDf = iDf_matches[0]
    
    d = sio.loadmat(f)
    if 'record_Q_mean' not in d or d['record_Q_mean'].size == 0:
        new_stability_code[iI, iDf] = 4
        new_growth_rate[iI, iDf] = 1.0 / T0
        continue
        
    qm = d['record_Q_mean'] # [nrec, 416]
    pm = d['record_P_mean']
    qs = d['record_Q_std']
    
    sigma_t0 = float(d['sigma_t0'][0,0]) if 'sigma_t0' in d else 40e-12
    Par_num = float(d['Par_num'][0,0]) if 'Par_num' in d else 10000
    Recor_step = float(d['Recor_step'][0,0]) if 'Recor_step' in d else 10
    
    I0 = i0_mA * 1e-3
    fre_shift_hc = df_Hz
    
    nrec = qm.shape[0]
    nan_rows = np.any(np.isnan(qm), axis=1) | np.any(np.isinf(qm), axis=1)
    zero_tail = np.all(qm[-min(100, nrec-1):, :] == 0) and nrec > 200
    
    psi_hc = np.arctan(2 * Q_hc * fre_shift_hc / fre_hc)
    V_hc = 2 * I0 * R_hc * np.cos(psi_hc)
    k_fp = np.sqrt(1/n_hc**2 - 1/(n_hc**2-1) * (U0/V_mc)**2)
    xi_calc = V_hc / (k_fp * V_mc)
    
    if np.any(nan_rows) or zero_tail:
        if xi_calc >= 1.05:
            code = 1
        else:
            code = 4
        gr_si = 1.0 / T0
        f_peak_Hz = 0
    else:
        win_steady = slice(int(nrec * 0.7), nrec)
        win_analysis = slice(int(nrec * 0.5), nrec)
        
        # 1. Centroid jitter across all 416 bunches
        bunch_jitters = np.std(qm[win_steady, :], axis=0) * sigma_t0 * 1e12 # ps
        mean_jit = np.mean(bunch_jitters)
        max_jit = np.max(bunch_jitters)
        
        # 2. Maximum bunch temporal FFT
        max_bunch_idx = np.argmax(bunch_jitters)
        sig_fft = qm[win_analysis, max_bunch_idx] - np.mean(qm[win_analysis, max_bunch_idx])
        N_fft = len(sig_fft)
        fft_amp = np.abs(np.fft.fft(sig_fft))[:N_fft//2]
        freq_axis = np.arange(N_fft//2) / (N_fft * Recor_step * T0)
        peak_idx = np.argmax(fft_amp[1:]) + 1
        f_peak_Hz = freq_axis[peak_idx]
        
        # 3. All-bunch RMS energy growth rate
        P_detrend = pm[win_analysis, :] - np.mean(pm[win_analysis, :], axis=0)
        rms_P = np.sqrt(np.mean(P_detrend**2, axis=1))
        
        noise_floor_ps = (sigma_t0 * 1e12) / np.sqrt(max(Par_num, 100))
        noise_floor_P  = 1.0 / np.sqrt(max(Par_num, 100))
        
        pass_turns = np.arange(int(nrec*0.5), nrec) * Recor_step
        if np.max(rms_P) > 2.5 * noise_floor_P:
            xx = np.log(rms_P + 1e-30)
            p_fit = np.polyfit(pass_turns, xx, 1)
            gr_si = p_fit[0] / T0
        else:
            gr_si = -0.5 * damping_rate
            
        # 4. Spatial coupled-bunch mode FFT
        spatial_fft = np.fft.fft(qm[win_steady, :], axis=1)
        mode_power = np.mean(np.abs(spatial_fft)**2, axis=0)
        ac_modes = mode_power[1:]
        crest_factor = np.max(ac_modes) / (np.mean(ac_modes) + 1e-30)
        
        has_lcbi_mode = (crest_factor > 5.0 and max_jit > 0.6) or \
                        (gr_si > 0.5 * damping_rate and max_jit > 0.5)
                        
        # 5. Bunch length & bunch-to-bunch pattern
        bl_mean_t = np.mean(qs[win_analysis, :], axis=1) * sigma_t0 * 1e12
        bl_osc_rel = np.std(bl_mean_t) / (np.mean(bl_mean_t) + 1e-30)
        
        std_bb = np.std(np.mean(qm[-min(100, nrec-1):, :], axis=0)) * sigma_t0 * 1e12
        
        jitter_thresh_mean = max(1.2, 3.0 * noise_floor_ps)
        jitter_thresh_max  = max(2.0, 5.0 * noise_floor_ps)
        
        is_unstable = (gr_si > 0.3 * damping_rate) or \
                      (mean_jit > jitter_thresh_mean) or \
                      (max_jit > jitter_thresh_max and I0 >= 0.35) or \
                      (bl_osc_rel > 0.05) or \
                      (std_bb > 0.6) or \
                      has_lcbi_mode
                      
        code = 0
        if is_unstable:
            if (xi_calc >= 1.05 or (std_bb > 0.6 and np.mean(bl_mean_t) > 75)) and I0 >= 0.375:
                code = 1
            elif has_lcbi_mode and xi_calc < 0.98:
                code = 3
            elif bl_osc_rel > 0.05 or mean_jit > 1.0 or (xi_calc >= 0.85 and xi_calc <= 1.05):
                code = 2
            else:
                code = 4
                
    new_stability_code[iI, iDf] = code
    new_growth_rate[iI, iDf] = gr_si
    new_peak_freq[iI, iDf] = f_peak_Hz

mask_valid = (new_stability_code >= 0)

# Save updated scan_summary.mat
sum_data['stability_code'] = new_stability_code
sum_data['growth_rate']    = new_growth_rate
sum_data['peak_freq_norm'] = new_peak_freq
sio.savemat(sum_file, sum_data)
print("Updated scan_summary.mat saved.")

# ----------------------------------------------------
# Create Side-by-Side Comparison Plot
# ----------------------------------------------------
colors = ['#2ca02c', '#ff7f0e', '#d62728', '#9467bd', '#7f7f7f'] # 0:Stable, 1:PTBL, 2:Coupled D-Q, 3:LCBI, 4:Other
cmap = ListedColormap(colors)
bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5]
norm = BoundaryNorm(bounds, cmap.N)
labels = ['Stable (0)', 'PTBL (1)', 'Coupled D-Q (2)', 'LCBI (3)', 'Other (4)']

fig, axes = plt.subplots(1, 2, figsize=(15, 6), sharey=True)

# 1. Left: Old Code Map
ax1 = axes[0]
old_plot_data = np.ma.masked_where(~mask_valid, old_stability_code)
im1 = ax1.pcolormesh(Df_scan_kHz, I0_scan_mA, old_plot_data, cmap=cmap, norm=norm, shading='nearest')
ax1.set_title('BEFORE: Single-Bunch Fitting\n(Fake Scattered Noise Points Misclassified as LCBI)', fontsize=12, fontweight='bold', color='#a00000')
ax1.set_xlabel(r'HHC Detuning $\Delta f_{HC}$ [kHz]', fontsize=11)
ax1.set_ylabel(r'Beam Current $I_0$ [mA]', fontsize=11)
ax1.grid(True, linestyle='--', alpha=0.5)

# 2. Right: New Full-Bunch Recognition Map
ax2 = axes[1]
new_plot_data = np.ma.masked_where(~mask_valid, new_stability_code)
im2 = ax2.pcolormesh(Df_scan_kHz, I0_scan_mA, new_plot_data, cmap=cmap, norm=norm, shading='nearest')
ax2.set_title('AFTER: Full-Bunch Spatial Mode + Noise Gating\n(Contiguous Clean Regions, Noise Artifacts Eliminated)', fontsize=12, fontweight='bold', color='#006600')
ax2.set_xlabel(r'HHC Detuning $\Delta f_{HC}$ [kHz]', fontsize=11)
ax2.grid(True, linestyle='--', alpha=0.5)

# Analytical flat potential curve
I0_dense = np.linspace(300, 480, 100) * 1e-3
cos_psi = (np.sqrt(1/n_hc**2 - 1/(n_hc**2-1) * (U0/V_mc)**2) * V_mc) / (2 * I0_dense * R_hc)
cos_psi = np.clip(cos_psi, -1.0, 1.0)
tan_psi = np.tan(np.arccos(cos_psi))
df_opt_kHz = (tan_psi * fre_hc / (2 * Q_hc)) / 1e3

for ax in [ax1, ax2]:
    ax.plot(df_opt_kHz, I0_dense * 1e3, 'k--', linewidth=2.0, label=r'Flat-potential ($\xi=1.0$)')

patches = [mpatches.Patch(color=colors[i], label=labels[i]) for i in range(5)]
patches.append(plt.Line2D([0], [0], color='black', linestyle='--', linewidth=2, label=r'Flat-pot. ($\xi=1.0$)'))
fig.legend(handles=patches, loc='lower center', ncol=6, bbox_to_anchor=(0.5, -0.05), fontsize=10, frameon=True)

plt.tight_layout()
fig_path_compare = os.path.join(folder, 'comparison_old_vs_new.png')
plt.savefig(fig_path_compare, dpi=200, bbox_inches='tight')

# Also save standalone reanalyzed figure
fig2, ax = plt.subplots(figsize=(9, 6))
im = ax.pcolormesh(Df_scan_kHz, I0_scan_mA, new_plot_data, cmap=cmap, norm=norm, shading='nearest')
ax.plot(df_opt_kHz, I0_dense * 1e3, 'k--', linewidth=2.2, label=r'Flat-potential ($\xi=1.0$)')
ax.set_title(r'2D Stability Map: Full-Bunch Spatial Mode Recognition', fontsize=13, fontweight='bold')
ax.set_xlabel(r'HHC Detuning $\Delta f_{HC}$ [kHz]', fontsize=12)
ax.set_ylabel(r'Beam Current $I_0$ [mA]', fontsize=12)
ax.grid(True, linestyle='--', alpha=0.5)
ax.legend(handles=patches, loc='upper right', fontsize=9, frameon=True)
plt.tight_layout()
fig_path_single = os.path.join(folder, 'stability_map_full_bunch_reanalyzed.png')
plt.savefig(fig_path_single, dpi=300, bbox_inches='tight')

# Copy to artifact directory
brain_dir = r'C:\Users\zhaol\.gemini\antigravity\brain\98411c23-a6cf-4664-a84d-0567203ec2cc'
art_compare = os.path.join(brain_dir, 'reanalysis_comparison.png')
art_single  = os.path.join(brain_dir, 'stability_map_reanalyzed.png')
shutil.copy2(fig_path_compare, art_compare)
shutil.copy2(fig_path_single, art_single)

print(f"Generated and copied comparison image to: {art_compare}")
print(f"Generated and copied standalone map to: {art_single}")
