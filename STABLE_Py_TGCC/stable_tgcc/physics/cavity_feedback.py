"""
Exact Python port of MATLAB STABLE cavity propagation and feedback controllers:
  - DFB_Set.m           -> dfb_set()
  - DFB_calc_properties -> dfb_calc_properties()
  - VoltageLoadCalc     -> voltage_load_calc()
  - MC_Cavity_Control.m -> mc_cavity_control()

Note on Direct RF Feedback (DFB) phase convention:
  `phase_shift = 0.0` [rad] represents ideal -180 deg negative feedback:
  `vg_drf = -DFB["gain"] * vc_delayed * exp(1j * DFB["phase_shift_app"])`
  where `phase_shift_app = psi_mc - phase_shift`.
"""

from __future__ import annotations
from typing import Dict, Any, Optional, Tuple
from collections import deque
import math
import cmath
import numpy as np
from scipy.signal import lfilter


def dfb_set(
    HALF: Dict[str, Any],
    gain: float = 0.2,
    phase_shift: float = 0.0,
    delay: int = 800,
    sample_num: int = 1,
    every: int = 1,
) -> Dict[str, Any]:
    """
    Initialize Direct RF Feedback (DFB) state dictionary (exact port of DFB_Set.m).
    `phase_shift = 0.0` [rad] is standard negative feedback.
    """
    every_i = max(1, int(round(every)))
    m_i = max(1, int(round(sample_num)))
    d_i = max(1, int(round(delay)))

    psi = float(HALF["det_angle_mc"])
    phase_shift_f = float(phase_shift)
    phase_shift_app = psi - phase_shift_f
    RL = float(HALF["R_mc"])
    QL = float(HALF["Q_mc"])
    RoverQ = RL / QL

    record_size = max(1, int(math.ceil(d_i / every_i)))
    vc_init = complex(HALF.get("Vrf_ideal", -HALF["V_mc"] * cmath.exp(-1j * (math.pi / 2.0 - HALF["fais_mc_whc"]))))

    # Negative feedback polarity built in so phase_shift = 0 is ideal negative feedback
    vg_drf_init = -float(gain) * vc_init * cmath.exp(1j * phase_shift_app)
    ig_dfb0 = vg_drf_init / (2.0 * RL * math.cos(psi)) * cmath.exp(-1j * psi)
    vg0_dfb = 2.0 * math.pi * RoverQ * ig_dfb0

    ig_mc_0 = complex(HALF.get("Ig_mc_0", 0j))
    ig_main0 = ig_mc_0 - ig_dfb0
    vg0_main = 2.0 * math.pi * RoverQ * ig_main0

    return {
        "enable": True,
        "gain": float(gain),
        "phase_shift": phase_shift_f,
        "d": d_i,
        "m": m_i,
        "every": every_i,
        "psi": psi,
        "phase_shift_app": phase_shift_app,
        "RL": RL,
        "QL": QL,
        "RoverQ": RoverQ,
        "record_size": record_size,
        "Ig_dfb0": ig_dfb0,
        "Vg0_dfb": vg0_dfb,
        "Ig_main0": ig_main0,
        "Vg0_main": vg0_main,
        "VcRecord": deque([vc_init] * record_size, maxlen=record_size),
        "Vrf_sample_buf": np.zeros(m_i, dtype=np.complex128),
        "dfbIndex": 0,
        "Ig_dfb": ig_dfb0,
        "Vg0": vg0_dfb,
    }


def m0_damper_set(
    HALF: Dict[str, Any],
    gain_kV_ps: float = 0.20,
    phase_deg: float = 0.0,
    delay_turns: int = 1,
    hp_alpha: float = 0.002,
    max_v_kV: float = 100.0,
) -> Dict[str, Any]:
    """
    Initialize MAX IV-style Mode-0 Longitudinal Damper Feedback on the Main Cavity.
    Measures turn-by-turn barycentric beam arrival time <tau> [ps], removes DC
    equilibrium offset via a 1-pole high-pass filter, computes the 90-deg
    synchrotron-quadrature velocity signal (with configurable phase shift), and
    injects a compensating generator current modulation dI_g into the Main Cavity.
    """
    psi = float(HALF["det_angle_mc"])
    RL = float(HALF["R_mc"])
    QL = float(HALF["Q_mc"])
    RoverQ = RL / QL
    T0 = float(HALF["T0"])

    # Estimate natural synchrotron phase advance per turn (omega_s * T0) for gain normalization
    v_mc = float(HALF["V_mc"])
    u0 = float(HALF["U0"])
    phi_s = float(HALF.get("fais_mc_whc", HALF.get("fais_nat", math.pi - math.asin(min(0.9999, u0 / max(v_mc, 1.0))))))
    qs_approx = math.sqrt(
        max(1e-10, (float(HALF["alpha_c"]) / (2.0 * math.pi)) * float(HALF["h"]) * v_mc * abs(math.cos(phi_s)) / float(HALF["E0"]))
    )
    omega_s_T0 = max(1e-5, 2.0 * math.pi * qs_approx)

    d_turns = max(1, int(round(delay_turns)))
    vg2ig_coef = cmath.exp(-1j * psi) / (2.0 * RL * max(1e-6, math.cos(psi)))
    roq_2pi = 2.0 * math.pi * RoverQ

    return {
        "enable": True,
        "gain_kV_ps": float(gain_kV_ps),
        "phase_rad": math.radians(float(phase_deg)),
        "delay_turns": d_turns,
        "hp_alpha": min(0.5, max(1e-5, float(hp_alpha))),
        "max_v_V": max(1.0, float(max_v_kV)) * 1e3,
        "omega_s_T0": omega_s_T0,
        "vg2ig_coef": vg2ig_coef,
        "roq_2pi": roq_2pi,
        "T0": T0,
        "dc_tau_ps": 0.0,
        "ac_prev_ps": 0.0,
        "initialized": False,
        "delay_buf": deque([0.0] * d_turns, maxlen=d_turns),
        "vg0_m0": 0j,
        "last_dv_kV": 0.0,
    }


def m0_damper_step(
    M0: Optional[Dict[str, Any]],
    tau_bary_ps: float,
) -> Optional[Dict[str, Any]]:
    """
    Update Mode-0 Damper state by one turn given the barycentric arrival time `tau_bary_ps` [ps],
    and compute the delayed Main Cavity generator kick `M0['vg0_m0']`.
    """
    if M0 is None or not M0.get("enable", False):
        return M0

    if not M0["initialized"]:
        M0["dc_tau_ps"] = float(tau_bary_ps)
        M0["ac_prev_ps"] = 0.0
        M0["initialized"] = True
        M0["vg0_m0"] = 0j
        return M0

    alpha = M0["hp_alpha"]
    dc = (1.0 - alpha) * M0["dc_tau_ps"] + alpha * float(tau_bary_ps)
    M0["dc_tau_ps"] = dc
    ac = float(tau_bary_ps) - dc

    # Quadrature velocity signal normalized to [ps] amplitude at f_s
    v_quad = (ac - M0["ac_prev_ps"]) / M0["omega_s_T0"]
    M0["ac_prev_ps"] = ac

    ph = M0["phase_rad"]
    sig_ps = v_quad * math.cos(ph) + ac * math.sin(ph)
    dv_V = M0["gain_kV_ps"] * 1e3 * sig_ps
    max_v = M0["max_v_V"]
    if dv_V > max_v:
        dv_V = max_v
    elif dv_V < -max_v:
        dv_V = -max_v

    M0["delay_buf"].appendleft(dv_V)
    dv_cmd_V = M0["delay_buf"].pop()
    M0["last_dv_kV"] = dv_cmd_V * 1e-3

    # Convert requested in-phase Main Cavity voltage modulation dV into per-bucket generator kick
    ig_m0 = dv_cmd_V * M0["vg2ig_coef"]
    M0["vg0_m0"] = M0["roq_2pi"] * ig_m0
    return M0


def dfb_calc_properties(
    DFB: Dict[str, Any],
    HALF: Dict[str, Any],
    Vc_val: Optional[complex] = None,
    Ig_val: Optional[complex] = None,
    Ig_dfb_val: Optional[complex] = None,
) -> Dict[str, Any]:
    """
    Calculate physical and diagnostic properties of Direct RF Feedback
    (exact port of DFB_calc_properties.m).
    """
    if Vc_val is None:
        Vc_val = complex(HALF.get("Vrf_ideal", HALF["V_mc"]))
    if Ig_val is None:
        Ig_val = complex(HALF.get("Ig_mc_0", 1.0 + 0j))
    if Ig_dfb_val is None:
        if "Ig_dfb" in DFB:
            Ig_dfb_val = complex(DFB["Ig_dfb"])
        else:
            phase_app = DFB["psi"] - DFB["phase_shift"]
            vg_drf = -DFB["gain"] * Vc_val * cmath.exp(1j * phase_app)
            Ig_dfb_val = vg_drf / (2.0 * DFB["RL"] * math.cos(DFB["psi"])) * cmath.exp(-1j * DFB["psi"])

    psi_dfb = cmath.phase(Vc_val) - cmath.phase(Ig_val)
    ig_main_val = Ig_val - Ig_dfb_val
    fac = abs(Ig_dfb_val) / max(abs(ig_main_val), 1e-30)
    alpha_db = 20.0 * math.log10(max(fac, 1e-30))
    gamma = fac / (1.0 - fac) if fac < 1.0 else float("inf")
    rs_eff = float(HALF["R_mc"]) / (1.0 + gamma * math.cos(psi_dfb))

    phase_app = DFB["psi"] - DFB["phase_shift"]
    vg_drf = -DFB["gain"] * Vc_val * cmath.exp(1j * phase_app)
    vg_main = complex(HALF.get("Vg_mc_init", 0j)) - vg_drf

    phi_s = float(HALF.get("fais_mc_whc", HALF.get("fais_nat", 2.7)))
    if "Vg_mc_init" in HALF:
        vc_fixed = vg_main / (1.0 + DFB["gain"] * cmath.exp(1j * phase_app))
    else:
        vc_fixed = Vc_val

    th_fixed = cmath.phase(vc_fixed) - cmath.phase(Vc_val) + phi_s
    vg_sum = -abs(vc_fixed) * math.cos(th_fixed)
    omega_s = 0.0
    if vg_sum > 0.0:
        omega_s = math.sqrt(HALF["alpha_c"] * HALF["w_rf"] * vg_sum / (HALF["E0"] * HALF["T0"]))

    return {
        "psi": psi_dfb,
        "alpha": alpha_db,
        "gamma": gamma,
        "Rs": rs_eff,
        "Vg_drf": vg_drf,
        "Vg_main": vg_main,
        "fs": omega_s / (2.0 * math.pi),
    }


def voltage_load_calc(
    V_load_0: complex,
    V_load_cpu: np.ndarray,
    TbAng_coef: complex,
    pattern: np.ndarray,
) -> Tuple[np.ndarray, complex]:
    """
    Exact, vectorized C-speed implementation of VoltageLoadCalc_matlab.m using `scipy.signal.lfilter`.
    """
    n = len(pattern)
    # Build per-bucket beam injection array `u` where u[ii] = V_load_cpu[j] if pattern[ii]==1 else 0
    u = np.zeros(n, dtype=np.complex128)
    u[pattern == 1] = V_load_cpu

    # Recurrence in VoltageLoadCalc_matlab.m:
    #   state_before[ii] = state_after[ii-1] * TbAng_coef
    #   state_after[ii]  = state_before[ii] + u[ii]
    #   V_load[ii]       = state_before[ii] (if pattern[ii]==1) else state_after[ii] (which is identical since u[ii]==0!)
    # Thus V_load[ii] == state_before[ii] for ALL ii, and state_after[ii] = state_before[ii] + u[ii]!
    # Linear filter for state_after: y[k] - TbAng_coef * y[k-1] = u[k], with y[-1] = V_load_0
    zi = np.array([V_load_0 * TbAng_coef], dtype=np.complex128)
    state_after, _ = lfilter([1.0], [1.0, -TbAng_coef], u, zi=zi)
    V_load = state_after - u
    return V_load, complex(state_after[-1])


def mc_cavity_control(
    PI: Optional[Dict[str, Any]],
    DFB: Optional[Dict[str, Any]],
    Vrf_ideal: complex,
    Vg_mc_track_0: complex,
    V_mc_load_0: complex,
    V_load_cpu: np.ndarray,
    TbAng_coef_mc: complex,
    pattern: np.ndarray,
    M0_DAMPER: Optional[Dict[str, Any]] = None,
) -> Tuple[np.ndarray, np.ndarray, complex, np.ndarray, complex, Optional[Dict[str, Any]], Optional[Dict[str, Any]]]:
    """
    Exact Python port of MC_Cavity_Control.m with O(1) deque FIFO buffers,
    plus optional MAX IV-style Mode-0 Damper Feedback modulation (`M0_DAMPER`).
    """
    n = len(pattern)
    Vg_mc_track = np.zeros(n, dtype=np.complex128)
    Vc_mc = np.zeros(n, dtype=np.complex128)

    # 1. Propagate beam loading across all buckets in one vectorized call
    V_load, V_mc_load_0 = voltage_load_calc(V_mc_load_0, V_load_cpu, TbAng_coef_mc, pattern)

    use_PI = (PI is not None) and ("m" in PI) and bool(PI.get("enable", True))
    use_DFB = (DFB is not None) and ("d" in DFB) and bool(DFB.get("enable", True))
    vg0_m0 = complex(M0_DAMPER["vg0_m0"]) if (M0_DAMPER is not None and M0_DAMPER.get("enable", False)) else 0j

    if use_PI and use_DFB and PI["piIndex"] == 0 and "Vg0_main" in DFB:
        PI["Vg0"] = DFB["Vg0_main"]
        PI["Ig"] = DFB["Ig_main0"]

    # Fast path when neither PI nor DFB updates every bucket dynamically
    if (not use_PI) and (not use_DFB):
        if PI is not None and "Vg0" in PI:
            vg0_const = complex(PI["Vg0"])
        elif DFB is not None and "Vg0_main" in DFB:
            vg0_const = complex(DFB["Vg0_main"])
        else:
            vg0_const = 0j
        u_gen = np.full(n, vg0_const + vg0_m0, dtype=np.complex128)
        zi = np.array([Vg_mc_track_0 * TbAng_coef_mc], dtype=np.complex128)
        Vg_mc_track, _ = lfilter([1.0], [1.0, -TbAng_coef_mc], u_gen, zi=zi)
        Vg_mc_track_0 = complex(Vg_mc_track[-1])
        Vc_mc = Vg_mc_track + V_load
        return Vc_mc, Vg_mc_track, Vg_mc_track_0, V_load, V_mc_load_0, PI, DFB

    # Extract locals for tight loop speed
    if use_PI and PI is not None:
        pi_m = int(PI["m"])
        pi_d = int(PI["d"])
        pi_rl = float(PI["RL"])
        pi_kp = float(PI["KP"])
        pi_ki = float(PI["KI"])
        pi_roq_2pi = 2.0 * math.pi * float(PI["RoverQ"])
        pi_fbid: deque = PI["Ig_FBid"]
        pi_buf: np.ndarray = PI["Vrf_mc_track"]
        pi_base_ig = complex(DFB["Ig_main0"]) if (use_DFB and DFB is not None and "Ig_main0" in DFB) else complex(PI["Ig0"])

    if use_DFB and DFB is not None:
        dfb_d = int(DFB["d"])
        dfb_m = int(DFB["m"])
        dfb_every = int(DFB["every"])
        dfb_gain = float(DFB["gain"])
        dfb_phase_exp = cmath.exp(1j * float(DFB["phase_shift_app"]))
        dfb_vg2ig_coef = cmath.exp(-1j * float(DFB["psi"])) / (2.0 * float(DFB["RL"]) * math.cos(float(DFB["psi"])))
        dfb_roq_2pi = 2.0 * math.pi * float(DFB["RoverQ"])
        dfb_rec: deque = DFB["VcRecord"]
        dfb_buf: np.ndarray = DFB["Vrf_sample_buf"]
        dfb_add_main = (not use_PI) and (PI is None or "Vg0" not in PI)
        dfb_vg0_main = complex(DFB.get("Vg0_main", 0j))

    # Vectorized C-speed fast path when loop latency d >= n (holds for all physical LLRF loops, e.g. d=704..800 >= h=416).
    # Works for ANY pi_m (e.g. 80) and ANY dfb_every / dfb_m (e.g. 1 or 208) even when h=416 is not divisible by pi_m!
    can_vectorize = (
        (not use_PI or (pi_d >= n and PI["piIndex"] >= pi_m))
        and (not use_DFB or (dfb_d >= n and DFB["dfbIndex"] >= dfb_m and len(dfb_rec) > (n // dfb_every + 2)))
    )
    if can_vectorize:
        u_gen = np.zeros(n, dtype=np.complex128)

        # 1. Construct PI generator kick vector across all n buckets
        if use_PI and PI is not None:
            pi_idx0 = PI["piIndex"]
            cur_vg0_pi = complex(PI["Vg0"])
            first_trig = -1
            for ii in range(n):
                p_idx = pi_idx0 + ii + 1
                if p_idx > pi_d and ((p_idx - pi_d) % pi_m == 0):
                    first_trig = ii
                    break
            if first_trig < 0:
                u_gen += cur_vg0_pi
            else:
                if first_trig > 0:
                    u_gen[:first_trig] += cur_vg0_pi
                pos = first_trig
                while pos < n:
                    if len(pi_fbid) > 0:
                        cur_vg0_pi = pi_roq_2pi * pi_fbid.popleft()
                    nxt = min(n, pos + pi_m)
                    u_gen[pos:nxt] += cur_vg0_pi
                    pos = nxt
            PI["Vg0"] = cur_vg0_pi
        elif PI is not None and "Vg0" in PI:
            u_gen += complex(PI["Vg0"])
        elif (not use_DFB) and (DFB is not None) and ("Vg0_main" in DFB):
            u_gen += complex(DFB["Vg0_main"])

        # 2. Construct DFB generator kick vector across all n buckets
        if use_DFB and DFB is not None:
            dfb_idx0 = DFB["dfbIndex"]
            cur_vg0_dfb = complex(DFB["Vg0"])
            base_add = dfb_vg0_main if dfb_add_main else 0j
            if dfb_every == 1 and dfb_idx0 >= dfb_d:
                # Vectorized 1-bucket step reading last n elements of dfb_rec in reverse
                rec_list = list(dfb_rec)
                vc_del_arr = np.array(rec_list[-1 : -1 - n : -1], dtype=np.complex128)
                vg_drf_arr = (-dfb_gain * dfb_phase_exp) * vc_del_arr
                ig_dfb_arr = vg_drf_arr * dfb_vg2ig_coef
                DFB["Ig_dfb"] = complex(ig_dfb_arr[-1])
                DFB["Vg0"] = complex(dfb_roq_2pi * ig_dfb_arr[-1])
                u_gen += dfb_roq_2pi * ig_dfb_arr + base_add
            else:
                first_dtrig = -1
                for ii in range(n):
                    d_idx = dfb_idx0 + ii + 1
                    if d_idx > dfb_d and ((d_idx - dfb_d) % dfb_every == 0):
                        first_dtrig = ii
                        break
                if first_dtrig < 0:
                    u_gen += cur_vg0_dfb + base_add
                else:
                    if first_dtrig > 0:
                        u_gen[:first_dtrig] += cur_vg0_dfb + base_add
                    pos = first_dtrig
                    while pos < n:
                        n_pushed_before = ((dfb_idx0 + pos) // dfb_every) - (dfb_idx0 // dfb_every)
                        vc_delayed = dfb_rec[-1 - n_pushed_before]
                        vg_drf = -dfb_gain * vc_delayed * dfb_phase_exp
                        ig_dfb = vg_drf * dfb_vg2ig_coef
                        DFB["Ig_dfb"] = ig_dfb
                        cur_vg0_dfb = dfb_roq_2pi * ig_dfb
                        nxt = min(n, pos + dfb_every)
                        u_gen[pos:nxt] += cur_vg0_dfb + base_add
                        pos = nxt
                DFB["Vg0"] = cur_vg0_dfb

        # 3. Single C-speed lfilter call for generator voltage across all n buckets
        if vg0_m0 != 0j:
            u_gen += vg0_m0
        zi = np.array([Vg_mc_track_0 * TbAng_coef_mc], dtype=np.complex128)
        Vg_mc_track, _ = lfilter([1.0], [1.0, -TbAng_coef_mc], u_gen, zi=zi)
        Vg_mc_track_0 = complex(Vg_mc_track[-1])
        Vc_mc = Vg_mc_track + V_load

        # 4. Update PI sampling & controller law
        if use_PI and PI is not None:
            pi_idx0 = PI["piIndex"]
            rem0 = pi_idx0 % pi_m
            # Reconstruct contiguous stream from the unconsumed prefix of pi_buf + Vc_mc
            stream_pi = np.concatenate((pi_buf[:rem0], Vc_mc)) if rem0 > 0 else Vc_mc
            n_full = len(stream_pi) // pi_m
            if n_full > 0:
                blocks = stream_pi[: n_full * pi_m].reshape(n_full, pi_m)
                means = np.mean(blocks, axis=1)
                for b in range(n_full):
                    di = (Vrf_ideal - complex(means[b])) / pi_rl
                    PI["Integral"] += pi_ki * di
                    PI["DIg"] = pi_kp * di + PI["Integral"]
                    PI["Ig"] = pi_base_ig + PI["DIg"]
                    pi_fbid.append(PI["Ig"])
                pi_buf[:] = blocks[-1]
            rem_end = len(stream_pi) % pi_m
            if rem_end > 0:
                pi_buf[:rem_end] = stream_pi[-rem_end:]
            PI["piIndex"] = pi_idx0 + n

        # 5. Update DFB sampling & FIFO buffer
        if use_DFB and DFB is not None:
            dfb_idx0 = DFB["dfbIndex"]
            if dfb_m == 1 and dfb_every == 1:
                for val in Vc_mc:
                    dfb_rec.appendleft(complex(val))
                dfb_buf[0] = Vc_mc[-1]
            elif dfb_m == dfb_every:
                rem0 = dfb_idx0 % dfb_every
                stream_dfb = np.concatenate((dfb_buf[:rem0], Vc_mc)) if rem0 > 0 else Vc_mc
                n_full = len(stream_dfb) // dfb_every
                if n_full > 0:
                    dblocks = stream_dfb[: n_full * dfb_every].reshape(n_full, dfb_every)
                    dmeans = np.mean(dblocks, axis=1)
                    for b in range(n_full):
                        dfb_rec.appendleft(complex(dmeans[b]))
                    dfb_buf[:] = dblocks[-1]
                rem_end = len(stream_dfb) % dfb_every
                if rem_end > 0:
                    dfb_buf[:rem_end] = stream_dfb[-rem_end:]
            else:
                for ii in range(n):
                    d_idx = dfb_idx0 + ii + 1
                    s_idx = (d_idx - 1) % dfb_m
                    dfb_buf[s_idx] = Vc_mc[ii]
                    if d_idx % dfb_every == 0:
                        dfb_rec.appendleft(complex(np.mean(dfb_buf)))
            DFB["dfbIndex"] = dfb_idx0 + n

        return Vc_mc, Vg_mc_track, Vg_mc_track_0, V_load, V_mc_load_0, PI, DFB

    for ii in range(n):
        # 1. PI trigger
        if use_PI and PI is not None:
            PI["piIndex"] += 1
            pi_idx = PI["piIndex"]
            if pi_idx > pi_d and ((pi_idx - pi_d) % pi_m == 0) and len(pi_fbid) > 0:
                PI["Vg0"] = pi_roq_2pi * pi_fbid.popleft()
            vg0_pi = PI["Vg0"]
        elif PI is not None and "Vg0" in PI:
            vg0_pi = PI["Vg0"]
        elif (not use_DFB) and (DFB is not None) and ("Vg0_main" in DFB):
            vg0_pi = DFB["Vg0_main"]
        else:
            vg0_pi = 0j

        # 2. DFB trigger
        if use_DFB and DFB is not None:
            DFB["dfbIndex"] += 1
            dfb_idx = DFB["dfbIndex"]
            if dfb_idx > dfb_d and ((dfb_idx - dfb_d) % dfb_every == 0):
                vc_delayed = dfb_rec[-1]
                vg_drf = -dfb_gain * vc_delayed * dfb_phase_exp
                ig_dfb = vg_drf * dfb_vg2ig_coef
                DFB["Ig_dfb"] = ig_dfb
                DFB["Vg0"] = dfb_roq_2pi * ig_dfb
            vg0_dfb = DFB["Vg0"] + (dfb_vg0_main if dfb_add_main else 0j)
        else:
            vg0_dfb = 0j

        # 3. Total generator voltage & cavity voltage at bucket ii
        Vg_mc_track_0 = Vg_mc_track_0 * TbAng_coef_mc + (vg0_pi + vg0_dfb + vg0_m0)
        Vg_mc_track[ii] = Vg_mc_track_0
        vc_cur = Vg_mc_track_0 + V_load[ii]
        Vc_mc[ii] = vc_cur

        # 4. Update PI sampling & controller law
        if use_PI and PI is not None:
            rem = PI["piIndex"] % pi_m
            if rem == 0:
                pi_buf[pi_m - 1] = vc_cur
                vrf_mc_mean = complex(np.mean(pi_buf))
                di = (Vrf_ideal - vrf_mc_mean) / pi_rl
                PI["Integral"] += pi_ki * di
                PI["DIg"] = pi_kp * di + PI["Integral"]
                PI["Ig"] = pi_base_ig + PI["DIg"]
                pi_fbid.append(PI["Ig"])
            else:
                pi_buf[rem - 1] = vc_cur

        # 5. Update DFB sampling & FIFO buffer
        if use_DFB and DFB is not None:
            s_idx = (DFB["dfbIndex"] - 1) % dfb_m
            dfb_buf[s_idx] = vc_cur
            if DFB["dfbIndex"] % dfb_every == 0:
                valid_len = min(DFB["dfbIndex"], dfb_m)
                mean_vc = complex(np.mean(dfb_buf[:valid_len]))
                dfb_rec.appendleft(mean_vc)

    return Vc_mc, Vg_mc_track, Vg_mc_track_0, V_load, V_mc_load_0, PI, DFB
