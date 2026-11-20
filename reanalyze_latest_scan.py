import os, glob, re, sys
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.patches as mpatches

def reanalyze_folder(folder):
    print(f"Target folder: {folder}")
    files = glob.glob(os.path.join(folder, 'SOLEIL_*.mat'))
    if not files:
        print("No files found!")
        return
    print(f"Found {len(files)} raw data files.")

    # 1. Parse grid coordinates from filenames
    records = []
    for f in files:
        m = re.search(r'I0(\d+)mA.*detune(-?\d+\.?\d*[eE]?[-+]?\d*)', os.path.basename(f))
        if m:
            i0_mA = int(m.group(1))
            df_Hz = float(m.group(2))
            records.append((i0_mA, df_Hz, f))

    I0_scan_mA = np.array(sorted(list(set(r[0] for r in records))))
    Df_scan_kHz = np.array(sorted(list(set(r[1]/1e3 for r in records))))
    nI0 = len(I0_scan_mA)
    nDf = len(Df_scan_kHz)

    print(f"Grid: {nI0} currents ({I0_scan_mA[0]} - {I0_scan_mA[-1]} mA) x {nDf} detunings ({Df_scan_kHz[0]:.1f} - {Df_scan_kHz[-1]:.1f} kHz)")

    stability_code_raw = np.full((nI0, nDf), -1, dtype=int)
    stability_code_new = np.full((nI0, nDf), -1, dtype=int)
    growth_rate_mat    = np.full((nI0, nDf), np.nan)
    peak_freq_mat      = np.full((nI0, nDf), np.nan)

    # Machine constants for SOLEIL II
    C = 353.98
    cspeed = 299792458
    T0 = C / cspeed
    h = 416
    tau_z = 13.7e-3
    damping_rate = 1.0 / tau_z
    R_hc = 1.86e6
    Q_hc = 31e3
    n_hc = 4
    U0 = 470e3
    V_mc = 1.7e6
    fre_hc = (h / T0) * n_hc

    for idx, (i0_mA, df_Hz, filepath) in enumerate(records):
        iI = np.where(np.abs(I0_scan_mA - i0_mA) < 1e-2)[0][0]
        iDf = np.where(np.abs(Df_scan_kHz - df_Hz/1e3) < 1e-2)[0][0]

        try:
            d = sio.loadmat(filepath)
        except Exception as e:
            print(f"Error loading {filepath}: {e}")
            continue

        if 'stability_code_point' in d:
            stability_code_raw[iI, iDf] = int(d['stability_code_point'][0,0])

        if 'record_Q_mean' not in d or d['record_Q_mean'].size == 0:
            stability_code_new[iI, iDf] = 4
            growth_rate_mat[iI, iDf] = 1.0 / T0
            continue

        qm = d['record_Q_mean']
        pm = d['record_P_mean']
        qs = d['record_Q_std']
        nrec = qm.shape[0]

        sigma_t0 = float(d['sigma_t0'][0,0]) if 'sigma_t0' in d else 40e-12
        Par_num = float(d['Par_num'][0,0]) if 'Par_num' in d else 10000
        Recor_step = float(d['Recor_step'][0,0]) if 'Recor_step' in d else 10

        I0 = i0_mA * 1e-3
        fre_shift_hc = df_Hz

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

            bunch_jitters = np.std(qm[win_steady, :], axis=0) * sigma_t0 * 1e12
            centroid_jitter_ps = np.mean(bunch_jitters)
            max_bunch_jitter_ps = np.max(bunch_jitters)

            max_bunch_idx = np.argmax(bunch_jitters)
            sig_fft = qm[win_analysis, max_bunch_idx] - np.mean(qm[win_analysis, max_bunch_idx])
            N_fft = len(sig_fft)
            fft_amp = np.abs(np.fft.fft(sig_fft))[:N_fft//2]
            freq_axis = np.arange(N_fft//2) / (N_fft * Recor_step * T0)
            peak_idx = np.argmax(fft_amp[1:]) + 1
            f_peak_Hz = freq_axis[peak_idx]

            P_detrend = pm[win_analysis, :] - np.mean(pm[win_analysis, :], axis=0)
            rms_P = np.sqrt(np.mean(P_detrend**2, axis=1))

            noise_floor_ps = (sigma_t0 * 1e12) / np.sqrt(max(Par_num, 100))
            noise_floor_P  = 1.0 / np.sqrt(max(Par_num, 100))

            pass_turns = np.arange(int(nrec*0.5), nrec) * Recor_step
            if np.max(rms_P) > 2.5 * noise_floor_P:
                xx = np.log(rms_P + 1e-30)
                try:
                    p_fit = np.polyfit(pass_turns, xx, 1)
                    gr_si = p_fit[0] / T0
                except:
                    gr_si = 0.0
            else:
                gr_si = -0.5 * damping_rate

            spatial_fft = np.fft.fft(qm[win_steady, :], axis=1)
            mode_power = np.mean(np.abs(spatial_fft)**2, axis=0)
            ac_modes = mode_power[1:]
            crest_factor = np.max(ac_modes) / (np.mean(ac_modes) + 1e-30)

            has_HOM = True
            has_lcbi_mode = ((crest_factor > 5.0 and max_bunch_jitter_ps > 0.6) or
                             (gr_si > 0.5 * damping_rate and max_bunch_jitter_ps > 0.5))

            bl_mean_t = np.mean(qs[win_analysis, :], axis=1) * sigma_t0 * 1e12
            bl_osc_rel = np.std(bl_mean_t) / (np.mean(bl_mean_t) + 1e-30)

            std_bb = np.std(np.mean(qm[-min(100, nrec-1):, :], axis=0)) * sigma_t0 * 1e12

            jitter_thresh_mean = max(1.2, 3.0 * noise_floor_ps)
            jitter_thresh_max  = max(2.0, 5.0 * noise_floor_ps)

            is_unstable = ((gr_si > 0.3 * damping_rate) or
                           (centroid_jitter_ps > jitter_thresh_mean) or
                           (max_bunch_jitter_ps > jitter_thresh_max and I0 >= 0.35) or
                           (bl_osc_rel > 0.08) or
                           (std_bb > 0.8) or
                           has_lcbi_mode)

            code = 0
            if is_unstable:
                # Flat-potential Landau-damped island: bounded non-growing motion is STABLE
                if (0.90 <= xi_calc <= 1.04) and (gr_si <= 0.35 * damping_rate) and (centroid_jitter_ps < 4.5):
                    code = 0  # Stable island
                elif (xi_calc >= 1.05 or (std_bb > 0.8 and np.mean(bl_mean_t) > 75)) and I0 >= 0.375:
                    code = 1  # PTBL
                elif has_lcbi_mode and xi_calc < 0.98:
                    code = 3  # LCBI
                elif bl_osc_rel > 0.08 or centroid_jitter_ps > 1.5:
                    code = 2  # Coupled D-Q
                else:
                    code = 0  # Damped / stable

        stability_code_new[iI, iDf] = code
        growth_rate_mat[iI, iDf] = gr_si
        peak_freq_mat[iI, iDf] = f_peak_Hz

        if (idx + 1) % 50 == 0 or (idx + 1) == len(records):
            print(f"Processed {idx+1}/{len(records)} points...")

    # Save scan_summary.mat
    sum_data = {
        'Df_scan_kHz': Df_scan_kHz,
        'I0_scan_mA': I0_scan_mA,
        'stability_code': stability_code_new,
        'growth_rate': growth_rate_mat,
        'peak_freq_norm': peak_freq_mat
    }
    sum_file = os.path.join(folder, 'scan_summary.mat')
    sio.savemat(sum_file, sum_data)
    print(f"Saved: {sum_file}")

    # Plot
    colors = ['#2ca02c', '#ff7f0e', '#d62728', '#9467bd', '#7f7f7f']
    cmap = ListedColormap(colors)
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5]
    norm = BoundaryNorm(bounds, cmap.N)
    labels = ['Stable (0)', 'PTBL (1)', 'Coupled D-Q (2)', 'LCBI (3)', 'Other (4)']

    fig, axes = plt.subplots(1, 2, figsize=(16, 7), sharey=True)

    # 1. Left: Raw Recorded Code
    ax1 = axes[0]
    mask_valid_raw = (stability_code_raw >= 0)
    data_raw = np.ma.masked_where(~mask_valid_raw, stability_code_raw)
    im1 = ax1.pcolormesh(Df_scan_kHz, I0_scan_mA, data_raw, cmap=cmap, norm=norm, shading='nearest')
    ax1.set_title(f'Raw Scan Output (Before Re-analysis)\n{os.path.basename(folder)}', fontsize=11, fontweight='bold')
    ax1.set_xlabel(r'$\Delta f_{\mathrm{HC}}$ (kHz)', fontsize=12)
    ax1.set_ylabel(r'Beam Current $I_0$ (mA)', fontsize=12)
    ax1.grid(True, linestyle=':', alpha=0.6)

    # 2. Right: Re-analyzed Code
    ax2 = axes[1]
    mask_valid_new = (stability_code_new >= 0)
    data_new = np.ma.masked_where(~mask_valid_new, stability_code_new)
    im2 = ax2.pcolormesh(Df_scan_kHz, I0_scan_mA, data_new, cmap=cmap, norm=norm, shading='nearest')
    ax2.set_title('Re-analyzed Stability Map\n(Protected Flat-Potential Stable Island)', fontsize=11, fontweight='bold', color='#0066cc')
    ax2.set_xlabel(r'$\Delta f_{\mathrm{HC}}$ (kHz)', fontsize=12)
    ax2.grid(True, linestyle=':', alpha=0.6)

    # Add Xi contours
    df_grid, i_grid = np.meshgrid(Df_scan_kHz*1e3, I0_scan_mA*1e-3)
    psi_grid = np.arctan(2 * Q_hc * df_grid / fre_hc)
    vhc_grid = 2 * i_grid * R_hc * np.cos(psi_grid)
    k_fp = np.sqrt(1/n_hc**2 - 1/(n_hc**2-1) * (U0/V_mc)**2)
    xi_grid = vhc_grid / (k_fp * V_mc)

    for ax in [ax1, ax2]:
        cs = ax.contour(Df_scan_kHz, I0_scan_mA, xi_grid, levels=[0.8, 0.9, 1.0, 1.05, 1.1], colors='white', linewidths=1.2, linestyles='--')
        ax.clabel(cs, fmt=r'$\xi=%.2f$', fontsize=9, inline=True)

    patches = [mpatches.Patch(color=colors[i], label=labels[i]) for i in range(len(labels))]
    fig.legend(handles=patches, loc='lower center', ncol=5, frameon=True, fontsize=11, bbox_to_anchor=(0.5, -0.02))

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.1)

    out_img = os.path.join(folder, 'stability_map_reanalyzed.png')
    plt.savefig(out_img, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Plot saved: {out_img}")

    # Copy to artifacts directory
    artifact_dir = r"C:\Users\zhaol\.gemini\antigravity\brain\98411c23-a6cf-4664-a84d-0567203ec2cc"
    dest_img = os.path.join(artifact_dir, "two_day_scan_reanalyzed.png")
    import shutil
    shutil.copyfile(out_img, dest_img)
    print(f"Copied to artifact: {dest_img}")

if __name__ == '__main__':
    latest_folder = r'scan_data\Scan_Df50-260kHz_I250-500mA_20260922_162237'
    reanalyze_folder(latest_folder)
