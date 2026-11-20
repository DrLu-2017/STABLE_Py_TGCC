import numpy as np
from machine_data import v2366_v3, load_TDR2_wf
from mbtrack2 import Beam, CavityResonator, LongitudinalMap, SynchrotronRadiation
from mbtrack2 import BeamMonitor, CavityMonitor, BunchSpectrumMonitor, BeamSpectrumMonitor
from mbtrack2 import WakePotential, WakePotentialMonitor, TunerLoop, ProportionalIntegralLoop, DirectFeedback, ProportionalLoop
from pathlib import Path

#%% Simulation Settings

# Flags
IDs = "open" # "open" or "close" or "close_phase2"
wake = True
P_FB = False
PI_FB = False
Tuner_FB = False
Direct_FB = True
spec = True
mpi = True
sweep_bool = False # sweep bunch/beam
remove_higher_mcf = False

## Filling
fill_type = "uniform" # "uniform" or "32b" or "single" or "gaps"
n_gaps = 0 # number of equidistant gap
gap_len = 0 # number of empty buckets per gap

# Calculation type
ContinuousTuning = False
Det_calc_Xi = False # If True, HC_det_end and HC_det are computed from xi_start and xi_end
Injection = False

# MC Phase calculation
HC_offset = True # calculcation using input parameters
HC_offset_Phasor = False # calculcation using simulation parameters
HC_offset_Phasor_power = False # calculcation using simulation parameters and power formula
impedance_offset = False # use bunch charge of bunch n°0 !
estimated_bunch_length = 40e-12
MCphase_update_during_tracking = False # use mean bunch length ?

# MC settings
Ncav_MC = 4
Set_optimal_detune = True
Injection_optimal_detune = False # keep optimal detune during current ramp

# HC settings
m_HC = 4
RoQ_HC = 30
Q0_HC = 31e3
QL_HC = 31e3
Ncav_HC = 2

# HOM settings
HOM_factor = 1.0

# detuned HC setting
Use_fixed_detuned_HC = False
fixed_detuned_HC = 425e3

# Settings 
Vc = 1.7e6
Itot = 0.338 # Use for fixed detuning or CT
HC_det = 67.88e3 #Hz   # = 71.6° # Use for fixed detuning 
HC_det_end = 3e3
MC_det = -35e3
xi_start = 1.0
xi_end = 1.1
Itot_final = 0.5

if Det_calc_Xi:
    if isinstance(xi_start, (int,float)):
        iter_values = [xi_start]
    else:
        iter_values = xi_start
else:
    if isinstance(HC_det, (int,float)):
        iter_values = [HC_det]
    else:
        iter_values = HC_det
N_iter = len(iter_values)

# FB Settings
fb_gain = [0.01,1000] # fb Gain for IQ components of Vc
fb_sample_num = 208 # 2*2*2*2*3*3*3, mean Vc value of this number is used for Vg control
fb_every = 208 # 192ns (assumption),corresponding to process speed of the Feedback system
fb_delay = 704 #int(2e-6/ring.T1)
fb_IIR_cutoff = 0 # Cutoff frequency of the IIR filter in [Hz]. If 0, cutoff frequency is infinity.
directFB_gain = 0.2
directFB_phaseShift = 0/180*np.pi
tuner_gain = 0.01
PFB_gainA = 0.01
PFB_gainP = 0.01
PFB_delay = 1

# Sweep settings
f0 = 400 # Hz
f1 = 2000 # Hz
t1 = 5000 # turns
level = 1e3 # V
plane = "tau"
bunch_to_sweep = None # None = all bunch

# Simulation Settings
mp = int(1e5)
tot_turns = int(300e3)
delay_before = int(0e3)
delay_after = int(0e3)
spec_step = int(50e3)
delay_before_sweep = int(50e3)
n_bin_wake = 80
tot_turns_iter = int(tot_turns*N_iter)

# Path
path = Path(__file__)

#%% Simulation setup

ring = v2366_v3(IDs=IDs, V_RF=Vc, HC_power=0)

if remove_higher_mcf:
    ring.mcf_order = ring.mcf_order[-1]

long = LongitudinalMap(ring)
rad = SynchrotronRadiation(ring)

#%% Setup monitors

spec_num = int(tot_turns_iter/spec_step)

name = path.parts[-1][:-3]

MCmon = CavityMonitor("MC", ring, file_name=name, save_every=100,
                 buffer_size=100, total_size=tot_turns_iter/100, mpi_mode=mpi)

HCmon = CavityMonitor("HHC", ring, file_name=None, save_every=100,
                 buffer_size=100, total_size=tot_turns_iter/100, mpi_mode=mpi)

bbmon = BeamMonitor(h=ring.h, file_name=None, save_every=10, buffer_size=100,
                 total_size=tot_turns_iter/10, mpi_mode=mpi)

if spec:
    beamspec = BeamSpectrumMonitor(ring, save_every=spec_step, buffer_size=1, total_size=spec_num, dim="tau",
                     n_fft=None, file_name=None, mpi_mode=mpi)

    bunchspec = BunchSpectrumMonitor(ring, bunch_number=0, mp_number=mp, sample_size=1e3, save_every=spec_step,
                     buffer_size=1, total_size=spec_num, dim="tau", n_fft=None,
                     file_name=None, mpi_mode=mpi, higher_orders=True)

if wake:
    wakemon = WakePotentialMonitor(0, "Wlong", n_bin_wake, save_every=100, buffer_size=100,
                                   total_size=int(tot_turns_iter/100), mpi_mode=mpi)

if Use_fixed_detuned_HC:
    HC_fixed_mon = CavityMonitor("HC_fixed", ring, file_name=None, save_every=100,
                 buffer_size=100, total_size=tot_turns_iter/100, mpi_mode=mpi)

for jj, val in enumerate(iter_values):
    
    if Det_calc_Xi:
        xi_start = val
    else:
        HC_det = val

    #%% Cavity setup
    
    m = 1
    Rs = 5e6
    Q = 35.7e3
    QL = 6e3
    detune = MC_det
    Ncav = Ncav_MC
    MC = CavityResonator(ring, m, Rs, Q, QL, detune, Ncav=Ncav)
    
    m = m_HC
    Rs = RoQ_HC * Q0_HC
    Q = Q0_HC
    QL = Q0_HC
    Ncav = Ncav_HC
    HHC = CavityResonator(ring, m, Rs, Q, QL, detune, Ncav=Ncav)
    HHC.Vg = 0
    HHC.theta_g = 0
    HHC.detune = HC_det
    
    #%% HOM 
    HOM = np.load("/ccc/work/cont003/soleil/gamelina/Data_mbtrack2/SOLEIL_II_HOM_MC_fit.npz")
    HOM_list = []
    nHOM = len(HOM["saved_Rs"])
    for i in range(nHOM):
        fr = (HOM["saved_M"][i]*ring.h + HOM["saved_mu"][i] + ring.synchrotron_tune(Vc)) * ring.f0
        m = fr/ring.f1
        HOMt = CavityResonator(ring, m, HOM["saved_Rs"][i]*HOM_factor, HOM["saved_Q"][i], HOM["saved_Q"][i], 100, Ncav=1)
        HOMt.Vg = 0
        HOMt.theta_g = 0
        HOMt.fr = fr
        HOM_list.append(HOMt)

    if Use_fixed_detuned_HC:
        HC_fixed = CavityResonator(ring, m, Rs, Q, QL, fixed_detuned_HC, Ncav=1)
        HC_fixed.Vg = 0
        HC_fixed.theta_g = 0
    
    #%% Wakes
    
    if wake:
        if IDs == "open":
            ID_wake = "open"
        else:
            ID_wake = "close"
        wf = load_TDR2_wf(f"TDR2.1_ID{ID_wake}")
        wf.drop(["Zxdip","Wxdip","Zydip","Wydip"])
        wp = WakePotential(ring, wf, n_bin=n_bin_wake)
    
    #%% Beam setup
    
    bb = Beam(ring)
    if fill_type == "uniform":
        filling = np.ones(ring.h)*Itot/ring.h
        nb = ring.h
    elif fill_type == "32b":
        filling = np.zeros(ring.h)
        filling[0::13] = Itot/32
        nb = 32
    elif fill_type == "single":
        filling = np.zeros(ring.h)
        filling[0] = Itot
        nb = 1
    elif fill_type == "gaps":
        filling = np.ones(ring.h)
        gap_var = int(ring.h/n_gaps)
        for i in range(n_gaps):
            filling[gap_var*(i+1) - gap_len:gap_var*(i+1)] = 0
        nb = ring.h - n_gaps*gap_len
        filling *= Itot/nb
    bb.init_beam(filling, mp_per_bunch=mp, track_alive=False, mpi=mpi)
    
    current_step = (Itot_final-Itot)/(tot_turns-delay_after)/nb
    
    #%% RF setup
    
    MC.Vc = Vc
    
    if Det_calc_Xi:
        HC_det = Itot*HHC.Rs/HHC.Q*ring.f1/Vc*HHC.m**2/xi_start
        HHC.detune = HC_det
        HC_det_end = Itot*HHC.Rs/HHC.Q*ring.f1/Vc*HHC.m**2/xi_end
    
    delta = 0
    if HC_offset or HC_offset_Phasor or HC_offset_Phasor_power:
        delta += HHC.Vb(Itot)*np.cos(HHC.psi)
    if impedance_offset:
        delta += bb[0].charge * wf.Wlong.loss_factor(estimated_bunch_length)
    MC.theta = np.arccos((ring.U0 + delta)/MC.Vc)
    
    if Set_optimal_detune:
        MC.set_optimal_detune(Itot)
    MC.set_generator(Itot)
    
    if P_FB:
        PFB = ProportionalLoop(ring, MC, PFB_gainA, PFB_gainP, PFB_delay)
        MC.feedback.append(PFB)
    if PI_FB and not Direct_FB:
        PI = ProportionalIntegralLoop(ring, MC, fb_gain, fb_sample_num, fb_every, fb_delay, IIR_cutoff=fb_IIR_cutoff)
        MC.feedback.append(PI)
    if Tuner_FB:
        Tuner = TunerLoop(ring, MC, tuner_gain)
        MC.feedback.append(Tuner)
    if Direct_FB:
        DFB = DirectFeedback(ring=ring, cav_res=MC, gain=fb_gain, sample_num=fb_sample_num,
                         every=fb_every, delay=fb_delay, IIR_cutoff=fb_IIR_cutoff,
                         DFB_gain=directFB_gain, DFB_phase_shift=directFB_phaseShift)
        MC.feedback.append(DFB)
    
    MC.init_phasor(bb)
    HHC.init_phasor(bb)
    if Use_fixed_detuned_HC:
        HC_fixed.init_phasor(bb)
    
    det_step = np.abs(HC_det - HC_det_end)/(tot_turns-delay_before-delay_after)
    
    #%% Sweep
    if sweep_bool:
        from mbtrack2 import Sweep
        sweep = Sweep(ring, f0, f1, ring.T0*t1, level, plane, bunch_to_sweep)
    
    #%% Tracking start
    
    for i in range(int(tot_turns+1)):
        if (i % 1000 == 0):
            if mpi:
                if ((bb.mpi.rank == 0) is True):
                    print(i)
            else:
                print(i)
    
        long.track(bb)
        rad.track(bb)
    
        if wake:
            wp.track(bb)
    
        if mpi:
            bb.mpi.share_distributions(bb)
    
        MC.track(bb)
        if i > delay_before:
            HHC.track(bb)
    
        if i > 50e3:
            for HOM in HOM_list:
                HOM.track(bb)

        if sweep_bool and (i > delay_before_sweep):
            sweep.track(bb)
    
        if Use_fixed_detuned_HC:
            HC_fixed.track(bb)
    
        if ContinuousTuning:
            if i < (tot_turns-delay_after):
                HHC.detune -= det_step
    
        if Injection and (i < (tot_turns-delay_after)):
            for bunch in bb:
                bunch.current += current_step
            if Injection_optimal_detune:
                MC.set_optimal_detune(bb.current)
            MC.set_generator(bb.current)
    
        if MCphase_update_during_tracking:
            delta = 0
            if HC_offset:
                delta += HHC.Vb(Itot)*np.cos(HHC.psi)
            elif HC_offset_Phasor:
                Vc = np.abs(HHC.cavity_phasor)
                theta = np.angle(HHC.cavity_phasor)
                delta += Vc*np.cos(theta)
            elif HC_offset_Phasor_power:
                Vc = np.abs(HHC.cavity_phasor)
                delta += Vc**2 / (2*HHC.Rs*bb.current)
    
            MC.theta = np.arccos((ring.U0 + delta)/MC.Vc)
    
        bbmon.track(bb)
        MCmon.track(bb, MC)
        HCmon.track(bb, HHC)
        if spec:
            bunchspec.track(bb)
            beamspec.track(bb)
        if wake:
            wakemon.track(bb, wp)
        if Use_fixed_detuned_HC:
            HC_fixed_mon.track(bb, HC_fixed)

bbmon.close()

