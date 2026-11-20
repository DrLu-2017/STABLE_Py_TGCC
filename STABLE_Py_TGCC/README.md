# STABLE_Py_TGCC — Multi-A100 GPU Server-Client RF & Longitudinal Beam Dynamics Suite

`STABLE_Py_TGCC` is a standalone Python Client-Server port of the verified MATLAB `STABLE` (`STABLE_App.m`) application, engineered for multi-GPU deployment on **TGCC Irene NVIDIA A100** nodes while keeping the original MATLAB codebase completely untouched.

---

## 1. Key Capabilities

1. **100% Verified Physics Parity with MATLAB `STABLE`**:
   - Self-consistent Fokker-Planck potential-well bunch form factor iteration $F_n(\sigma_\tau)$ (`calc_bunch_form_factor_eq`).
   - Exact Robinson reactive beam-loading compensation detuning (`compute_optimal_mc_detune`).
   - Independent PI and Direct RF Feedback (`DFB`) cavity control with standard negative-feedback polarity (`phase_shift = 0°` represents $-180^\circ$ negative feedback).
   - Multi-bunch intra-bunch short-range wake convolution + multi-HOM coupled-bunch tracking.
   - T. He et al. (PRAB 2024) analytical D-mode Robinson instability pre-screening & LLRF Bode/Nyquist stability margin analyzer.

2. **Multi-A100 GPU Backend (`CuPy` -> `PyTorch CUDA` -> `NumPy CPU`)**:
   - Automatically detects all available NVIDIA A100 GPUs on a TGCC node and dispatches 2D `(I0, Df_HC)` scan points in parallel across GPUs.

3. **Encrypted Multi-User Management (`SQLite` + `PBKDF2-HMAC-SHA256` + `AES/Fernet`)**:
   - Each user signs in with a master password to unlock their personal encrypted TGCC SSH credentials (`ssh_host`, `ssh_username`, encrypted password/passphrase, `-A` project account, `-q a100` partition) and personal `MachineConfig` / `ScanConfig` presets.

4. **Hybrid TGCC Server-Client Communication**:
   - **Local Server Mode**: Runs an embedded background RPC server on `127.0.0.1` for local testing and plotting.
   - **Remote TGCC SSH Mode**: Connects to `irene-fr.ccc.cea.fr` via SSH (supporting password, SSH key, and session 2FA/OTP token), deploys the server package via SFTP with one click, and streams live turn-by-turn waveforms and 2D stability map points over an SSH tunnel (`daemon_stream`) or via `ccc_msub` A100 batch jobs (`ccc_msub_batch`).

---

## 2. Quick Start

### Launch the Desktop GUI Client
```bash
cd d:\Lu\Projects\STABLE\STABLE_Py_TGCC
python run_client.py
```

### Run the Automated Verification Suite
```bash
cd d:\Lu\Projects\STABLE\STABLE_Py_TGCC
python -m unittest discover -s tests -p "test_*.py" -v
```
