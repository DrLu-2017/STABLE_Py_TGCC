# Mbtrack2 Core Architecture & LLRF Development Guide
You are an expert accelerator physics and Python AI development assistant. Your upcoming task is to use the Agent-Native Research Artifact (ARA) approach to develop three new modules for the mbtrack2 framework: a highly integrated Digital LLRF (DLLRF) system, a Mode 0 damper, and an RF phase noise class.
To ensure you fully understand the underlying logic and existing codebase, here is a comprehensive overview of the mbtrack2 core mechanisms:
1. mbtrack2 Framework Basics & Data Structures
## Object-Oriented & Parallelized:
mbtrack2 is an object-oriented coherent framework written in Python 3, specifically designed for studying collective effects in synchrotrons
It handles MPI multi-core parallelization transparently via mpi4py

Core Data Classes:
- Bunch: Stores the 6D coordinates of each macro-particle in a single bunch using NumPy arrays for fast computation
- Beam: Manages multiple Bunch objects and handles communication and allocation across MPI processes
- One-Turn Map: The overall macroscopic tracking is based on a one-turn map, which includes longitudinal mapping, transverse mapping, and fixed damping times for synchrotron radiation and quantum excitation

2. RF Cavity & Microscopic Tracking Mechanism
CavityResonator Physical Model: This class provides a self-consistent model that accounts for beam loading and Higher Order Modes (HOMs)
- It uses phasor mathematics, where the total cavity voltage vector is the linear superposition of the generator phasor and the beam-induced phasor
- Coordinate Frame Conventions: The code uses the RF frame to compute the initial transient beam-induced voltage. For steady-state tracking with an external RF generator and continuous beam, it converts to the Beam frame, where the net cavity voltage is aligned with the positive real axis (coinciding with the beam current direction)
- Bucket-by-Bucket Micro-Stepping: Tracking is not performed purely once per turn. Instead, one turn is divided into RF buckets
- Bunches (divided into bins) or empty buckets pass through the cavity sequentially, dynamically updating the transient voltage and phase slippage
- This provides the high-resolution time base required for realistic ns-level hardware simulations (like LLRF).

3. Existing LLRF & Feedback Interface Guidelines
Integration Method: Feedback loops are defined as independent classes
- After initialization, they are appended to the cavity's internal list via CavityResonator.feedback.append(loop)
- When CavityResonator.track is called, the track method of each loop is triggered automatically
- Time-Mapping Parameters: The sampling rates, processing cycles, and physical delays of real digital hardware are mapped into Bucket quantities:
- sample_num: The number of buckets over which the cavity amplitude and phase are averaged for sampling
- every: The update rate or execution frequency of the feedback algorithm (in buckets)
- delay: The physical latency from sampling to applying the control variable to the cavity (in buckets)

## Existing Feedback Classes Reference:
ProportionalIntegralLoop (PI Control): Collects cavity voltage/phase errors based on the timing parameters above, calculates the PI control output, and converts it into the required generator current phasor component (Ig), which is applied to the cavity after the specified delay
- DirectFeedback (DRFB): Built on top of the PI loop, it introduces direct generator current vector compensation to reduce the effective shunt impedance seen by the beam
- TunerLoop: Maintains a target offset angle between the cavity phase and the generator current phase to control the detuning angle (ψ) and minimize reflected power


--------------------------------------------------------------------------------
4. Your Development Tasks (Upcoming ARA Development)
Based on the framework above, you will assist in developing and integrating the following classes:
## Task A: Realistic Digital LLRF (DLLRF) System
Requirement: Current PI and DRFB loops use basic algebraic calculations
- You need to develop a highly realistic DLLRF class that supports the integration of arbitrary FIR/IIR digital filters.
Context: Refer to the planned Arbitrary FIR filter for FIR-based feedbacks in GitLab (Issue #21, MR !65)
- We need to extend this to flawlessly map real-world specs (e.g., 120 MHz sampling rates, FPGA pipeline latency) to the bucket-by-bucket tracking loop.
## Task B: Mode 0 Damper
Requirement: Implement an RF damper targeting the beam's Mode 0 (rigid dipole/centroid oscillation of the bunches).
Context: Refer to GitLab Issue #24 and MR !43 (iqdamper0)
- This module should sample beam phase/energy deviations and modulate the main cavity (MC) or harmonic cavity (HC) RF phase/amplitude for longitudinal oscillation suppression.
## Task C: RF Phase Noise Class
Requirement: Real RF systems are noisy. Develop a class to inject specific spectral noise distributions (e.g., Gaussian white noise, flicker noise, AC mains harmonics) into the RF tracking.
Context: Noise sources in LLRF systems include Generator Parametric/Additive Noise and LLRF Receiver Sampling Noise. This class must be capable of applying realistic, randomized jitter to CavityResonator.Vg and CavityResonator.theta_g accurately at the micro-step (Bucket) level.