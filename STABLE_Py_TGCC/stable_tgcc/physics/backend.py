"""
Unified GPU/CPU Array Backend for STABLE_Py_TGCC.
Automatically selects:
  1. CuPy (NVIDIA A100 CUDA backend)
  2. PyTorch CUDA (if CuPy is not installed on TGCC container)
  3. NumPy vectorized CPU fallback (for local development / testing without GPU)
"""

from __future__ import annotations
from typing import List, Dict, Any, Tuple
import numpy as np

try:
    import cupy as cp
    _HAS_CUPY = cp.cuda.is_available()
except Exception:
    cp = None
    _HAS_CUPY = False

try:
    import torch
    _HAS_TORCH_CUDA = torch.cuda.is_available()
except Exception:
    torch = None
    _HAS_TORCH_CUDA = False

from stable_tgcc.physics.cuda_kernels import get_compiled_cuda_lib, detect_nvcc_cuda_gpus


def detect_available_gpus() -> List[Dict[str, Any]]:
    """Return list of available NVIDIA GPUs (e.g., A100-SXM4-80GB on TGCC Topaze)."""
    gpus: List[Dict[str, Any]] = []
    if _HAS_CUPY and cp is not None:
        n_dev = cp.cuda.runtime.getDeviceCount()
        for i in range(n_dev):
            props = cp.cuda.runtime.getDeviceProperties(i)
            name = props["name"].decode("utf-8") if isinstance(props["name"], bytes) else str(props["name"])
            mem_gb = props["totalGlobalMem"] / (1024 ** 3)
            gpus.append({"id": i, "name": name, "memory_gb": round(mem_gb, 1), "backend": "cupy"})
        return gpus
    if _HAS_TORCH_CUDA and torch is not None:
        n_dev = torch.cuda.device_count()
        for i in range(n_dev):
            props = torch.cuda.get_device_properties(i)
            mem_gb = props.total_memory / (1024 ** 3)
            gpus.append({"id": i, "name": props.name, "memory_gb": round(mem_gb, 1), "backend": "pytorch_cuda"})
        return gpus
    nvcc_gpus = detect_nvcc_cuda_gpus()
    if nvcc_gpus:
        return nvcc_gpus
    return [{"id": 0, "name": "CPU (NumPy Vectorized)", "memory_gb": 0.0, "backend": "numpy"}]


class ArrayBackend:
    """
    Device-pinned array backend for running a single macro-particle tracking task
    on a specific GPU (`device_id`) or CPU.
    """

    def __init__(self, device_id: int = 0, prefer_cpu: bool = False):
        self.device_id = int(device_id)
        if not prefer_cpu and _HAS_CUPY and cp is not None:
            self.mode = "cupy"
            self.dev = cp.cuda.Device(self.device_id)
            self.dev.use()
            props = cp.cuda.runtime.getDeviceProperties(self.device_id)
            name = props["name"].decode("utf-8") if isinstance(props["name"], bytes) else str(props["name"])
            self.device_name = f"GPU:{self.device_id} ({name} [CuPy])"
        elif not prefer_cpu and _HAS_TORCH_CUDA and torch is not None:
            self.mode = "torch"
            self.torch_device = torch.device(f"cuda:{self.device_id}")
            self.device_name = f"GPU:{self.device_id} ({torch.cuda.get_device_name(self.device_id)} [PyTorch])"
        elif not prefer_cpu and get_compiled_cuda_lib() is not None:
            self.mode = "nvcc_cuda"
            self.cuda_lib = get_compiled_cuda_lib()
            nv_gpus = detect_nvcc_cuda_gpus()
            if self.device_id >= len(nv_gpus):
                self.device_id = 0
            gname = nv_gpus[self.device_id]["name"] if nv_gpus else "NVIDIA A100"
            self.device_name = f"GPU:{self.device_id} ({gname} [CUDA C++])"
        else:
            self.mode = "numpy"
            self.device_name = "CPU (NumPy Vectorized)"

    def to_device(self, arr: np.ndarray, dtype=np.float32):
        a = np.asarray(arr, dtype=dtype)
        if self.mode == "cupy":
            with self.dev:
                return cp.asarray(a)
        elif self.mode == "torch":
            return torch.as_tensor(a, device=self.torch_device)
        return a

    def to_numpy(self, arr) -> np.ndarray:
        if self.mode == "cupy":
            return cp.asnumpy(arr)
        elif self.mode == "torch":
            return arr.detach().cpu().numpy()
        return np.asarray(arr)

    def synchronize(self) -> None:
        if self.mode == "cupy":
            self.dev.synchronize()
        elif self.mode == "torch":
            torch.cuda.synchronize(self.torch_device)

    def randn(self, shape: Tuple[int, ...]):
        if self.mode == "cupy":
            with self.dev:
                return cp.random.standard_normal(shape, dtype=cp.float32)
        elif self.mode == "torch":
            return torch.randn(shape, device=self.torch_device, dtype=torch.float32)
        return np.random.standard_normal(shape).astype(np.float32)

    def max_abs(self, Q) -> float:
        if self.mode == "cupy":
            with self.dev:
                return float(cp.max(cp.abs(Q)).get())
        elif self.mode == "torch":
            return float(torch.max(torch.abs(Q)).item())
        return float(np.max(np.abs(Q)))

    def exp_and_col_sum(self, coef: complex, Q):
        """
        Compute exp_angle = exp(coef * Q) on device and return (exp_angle, sum_over_particles_numpy).
        Q has shape (Par_num, Bun_num).
        """
        c_re = np.float32(np.real(coef))
        c_im = np.float32(np.imag(coef))
        if self.mode == "cupy":
            with self.dev:
                arg = (c_re + 1j * c_im) * Q
                exp_a = cp.exp(arg)
                s_cpu = cp.asnumpy(cp.sum(exp_a, axis=0)).astype(np.complex128)
                return exp_a, s_cpu
        elif self.mode == "torch":
            mag = torch.exp(Q * float(c_re))
            ang = Q * float(c_im)
            exp_a = torch.complex(mag * torch.cos(ang), mag * torch.sin(ang))
            s_cpu = torch.sum(exp_a, dim=0).detach().cpu().numpy().astype(np.complex128)
            return exp_a, s_cpu
        else:
            arg = (c_re + 1j * c_im) * Q
            exp_a = np.exp(arg)
            s_cpu = np.sum(exp_a, axis=0).astype(np.complex128)
            return exp_a, s_cpu

    def apply_cavity_kick(self, P, v_load_filled_cpu: np.ndarray, exp_angle, imag_factor: float):
        """
        Apply cavity kick to normalized momentum P:
          V_kick = V_load_gpu / exp_angle
          P = P - real(V_kick) + imag(V_kick) * imag_factor
        """
        im_f = np.float32(imag_factor)
        if self.mode == "cupy":
            with self.dev:
                v_gpu = cp.asarray(v_load_filled_cpu.astype(np.complex64))[None, :]
                v_kick = v_gpu / exp_angle
                return P - cp.real(v_kick) + cp.imag(v_kick) * im_f
        elif self.mode == "torch":
            v_gpu = torch.as_tensor(v_load_filled_cpu.astype(np.complex64), device=self.torch_device).unsqueeze(0)
            v_kick = v_gpu / exp_angle
            return P - v_kick.real + v_kick.imag * float(im_f)
        else:
            v_kick = v_load_filled_cpu.astype(np.complex64)[None, :] / exp_angle
            return P - np.real(v_kick) + np.imag(v_kick) * im_f

    def short_range_wake_kick(
        self,
        Q,
        Dq: float,
        Wake_inter_cpu: np.ndarray,
        wake_kick_coef_filled: np.ndarray,
        ramp_factor: float,
    ) -> Tuple[bool, Any]:
        """
        Exact GPU/CPU implementation of BinNumCalZ.cu + conv2 short-range wake kick:
          Q_min = min(Q, axis=0)
          Q_new = round((Q - Q_min) / Dq)  (0-based bin index in [0, binnum-1])
          bin_num_q[b, j] = count of particles in bunch j with Q_new == b
          kick_conv = conv(bin_num_q[:, j], Wake_inter[0:binnum])[:binnum] * wake_kick_coef[j] * ramp_factor
          wake_kick[p, j] = kick_conv[Q_new[p, j], j]
        Returns (ok, wake_kick_on_device).
        """
        inv_dq = np.float32(1.0 / Dq)
        par_num, bun_num = Q.shape

        if self.mode == "cupy":
            with self.dev:
                q_min = cp.min(Q, axis=0, keepdims=True)
                q_idx = cp.rint((Q - q_min) * inv_dq).astype(cp.int32)
                max_idx = int(cp.max(q_idx).get())
                min_idx = int(cp.min(q_idx).get())
                binnum = max_idx + 1
                if min_idx < 0 or binnum < 1 or binnum > 2000 or binnum > len(Wake_inter_cpu):
                    return False, 0.0
                # Flattened 2D histogram via bincount: linear index = col * binnum + q_idx
                col_offsets = (cp.arange(bun_num, dtype=cp.int32) * binnum)[None, :]
                lin_idx = (q_idx + col_offsets).ravel()
                counts_flat = cp.bincount(lin_idx, minlength=bun_num * binnum).astype(cp.float32)
                # Shape (bun_num, binnum) -> transpose to (binnum, bun_num)
                bin_num_q = counts_flat.reshape(bun_num, binnum).T
                w_ker = cp.asarray(Wake_inter_cpu[:binnum].astype(np.float32))
                # 1D causal convolution along axis 0 via FFT
                n_fft = int(2 ** np.ceil(np.log2(2 * binnum - 1)))
                W_f = cp.fft.rfft(w_ker, n=n_fft)[:, None]
                B_f = cp.fft.rfft(bin_num_q, n=n_fft, axis=0)
                conv_full = cp.fft.irfft(B_f * W_f, n=n_fft, axis=0)[:binnum, :]
                scale = cp.asarray((wake_kick_coef_filled * ramp_factor).astype(np.float32))[None, :]
                kick_table = (conv_full * scale).T.ravel()  # (bun_num * binnum,)
                wake_kick = kick_table[lin_idx].reshape(par_num, bun_num)
                return True, wake_kick

        elif self.mode == "torch":
            q_min = torch.min(Q, dim=0, keepdim=True).values
            q_idx = torch.round((Q - q_min) * float(inv_dq)).to(torch.int64)
            max_idx = int(torch.max(q_idx).item())
            min_idx = int(torch.min(q_idx).item())
            binnum = max_idx + 1
            if min_idx < 0 or binnum < 1 or binnum > 2000 or binnum > len(Wake_inter_cpu):
                return False, 0.0
            col_offsets = (torch.arange(bun_num, device=self.torch_device, dtype=torch.int64) * binnum).unsqueeze(0)
            lin_idx = (q_idx + col_offsets).reshape(-1)
            counts_flat = torch.bincount(lin_idx, minlength=bun_num * binnum).to(torch.float32)
            bin_num_q = counts_flat.reshape(bun_num, binnum).T  # (binnum, bun_num)
            w_ker = torch.as_tensor(Wake_inter_cpu[:binnum].astype(np.float32), device=self.torch_device)
            n_fft = int(2 ** np.ceil(np.log2(2 * binnum - 1)))
            W_f = torch.fft.rfft(w_ker, n=n_fft).unsqueeze(1)
            B_f = torch.fft.rfft(bin_num_q, n=n_fft, dim=0)
            conv_full = torch.fft.irfft(B_f * W_f, n=n_fft, dim=0)[:binnum, :]
            scale = torch.as_tensor((wake_kick_coef_filled * ramp_factor).astype(np.float32), device=self.torch_device).unsqueeze(0)
            kick_table = (conv_full * scale).T.reshape(-1)
            wake_kick = kick_table[lin_idx].reshape(par_num, bun_num)
            return True, wake_kick

        else:
            q_min = np.min(Q, axis=0, keepdims=True)
            q_idx = np.rint((Q - q_min) * inv_dq).astype(np.int32)
            max_idx = int(np.max(q_idx))
            min_idx = int(np.min(q_idx))
            binnum = max_idx + 1
            if min_idx < 0 or binnum < 1 or binnum > 2000 or binnum > len(Wake_inter_cpu):
                return False, 0.0
            col_offsets = (np.arange(bun_num, dtype=np.int32) * binnum)[None, :]
            lin_idx = (q_idx + col_offsets).ravel()
            counts_flat = np.bincount(lin_idx, minlength=bun_num * binnum).astype(np.float32)
            bin_num_q = counts_flat.reshape(bun_num, binnum).T  # (binnum, bun_num)
            w_ker = Wake_inter_cpu[:binnum].astype(np.float32)
            n_fft = int(2 ** np.ceil(np.log2(2 * binnum - 1)))
            W_f = np.fft.rfft(w_ker, n=n_fft)[:, None]
            B_f = np.fft.rfft(bin_num_q, n=n_fft, axis=0)
            conv_full = np.fft.irfft(B_f * W_f, n=n_fft, axis=0)[:binnum, :].astype(np.float32)
            scale = (wake_kick_coef_filled * ramp_factor).astype(np.float32)[None, :]
            kick_table = (conv_full * scale).T.ravel()
            wake_kick = kick_table[lin_idx].reshape(par_num, bun_num)
            return True, wake_kick

    def mean_and_std_cols(self, Q, P) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Compute per-bunch mean and std of Q and P on device and return as NumPy arrays."""
        if self.mode == "cupy":
            with self.dev:
                qm = cp.asnumpy(cp.mean(Q, axis=0))
                qs = cp.asnumpy(cp.std(Q, axis=0))
                pm = cp.asnumpy(cp.mean(P, axis=0))
                ps = cp.asnumpy(cp.std(P, axis=0))
                return qm, qs, pm, ps
        elif self.mode == "torch":
            qm = torch.mean(Q, dim=0).detach().cpu().numpy()
            qs = torch.std(Q, dim=0, unbiased=False).detach().cpu().numpy()
            pm = torch.mean(P, dim=0).detach().cpu().numpy()
            ps = torch.std(P, dim=0, unbiased=False).detach().cpu().numpy()
            return qm, qs, pm, ps
        else:
            return (
                np.mean(Q, axis=0),
                np.std(Q, axis=0),
                np.mean(P, axis=0),
                np.std(P, axis=0),
            )
