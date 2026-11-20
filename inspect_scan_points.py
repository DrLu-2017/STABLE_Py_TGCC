import glob, re, scipy.io as sio, numpy as np

files = sorted(glob.glob('scan_data/Scan_Df50-250kHz_I300-500mA_20260921_105745/SOLEIL_I0300mA_*.mat'),
               key=lambda x: int(re.search(r'detune(\d+)', x).group(1)))

print(f"Inspecting {len(files)} files for 300mA:")
header = f"{'Df(kHz)':>8} | {'code':>4} | {'gr_si':>10} | {'f_peak(Hz)':>11} | {'mean_jit':>10} | {'max_jit':>9} | {'b1_jit':>9} | {'std_bb':>8} | {'bl_osc':>8}"
print(header)
print('-'*len(header))

for f in files[8:30]: # detuning ~90 kHz to 200 kHz
    m = re.search(r'detune(\d+)', f)
    df = int(m.group(1)) / 1e3
    d = sio.loadmat(f)
    code = int(d['stability_code_point'][0,0])
    gr = float(d['gr_si'][0,0])
    fp = float(d['f_peak_Hz'][0,0])
    
    sigma_t0 = float(d['sigma_t0'][0,0])
    qm = d['record_Q_mean'] # 3000 x 416
    pm = d['record_P_mean']
    qs = d['record_Q_std']
    
    nrec = qm.shape[0]
    win_steady = slice(int(nrec*0.8), nrec)
    
    bunch_jit = np.std(qm[win_steady, :], axis=0) * sigma_t0 * 1e12
    mean_jit = np.mean(bunch_jit)
    max_jit = np.max(bunch_jit)
    b1_jit = bunch_jit[0]
    
    std_bb = np.std(np.mean(qm[-100:, :], axis=0)) * sigma_t0 * 1e12
    bl_mean = np.mean(qs[win_steady, :], axis=1)
    bl_osc = np.std(bl_mean) / np.mean(bl_mean)
    
    print(f"{df:8.1f} | {code:4d} | {gr:10.2f} | {fp:11.1f} | {mean_jit:10.3f} | {max_jit:9.3f} | {b1_jit:9.3f} | {std_bb:8.3f} | {bl_osc:8.4f}")
