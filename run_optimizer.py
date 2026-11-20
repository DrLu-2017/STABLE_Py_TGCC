"""
CLI runner for RFOptimizer called from STABLE MATLAB GUI or terminal.
Usage:
    python run_optimizer.py <I0_mA> [V_mc] [U0] [n_hc] [R_hc] [Q_hc]
Output:
    <detuning_khz>,<xi>,<bunch_length_factor>,<status_name>
"""

import sys
from pathlib import Path

# Ensure workspace is in path
sys.path.insert(0, str(Path(__file__).resolve().parent))

from stable_core.types import MachineParameters
from stable_core.optimizer import RFOptimizer


def main() -> None:
    if len(sys.argv) < 2:
        print("ERROR: missing I0_mA argument")
        sys.exit(1)

    try:
        I0_mA = float(sys.argv[1])
        V_mc = float(sys.argv[2]) if len(sys.argv) > 2 else 1.7e6
        U0 = float(sys.argv[3]) if len(sys.argv) > 3 else 470e3
        n_hc = int(float(sys.argv[4])) if len(sys.argv) > 4 else 4
        R_hc = float(sys.argv[5]) if len(sys.argv) > 5 else 1.86e6
        Q_hc = float(sys.argv[6]) if len(sys.argv) > 6 else 31000.0

        machine = MachineParameters(
            I0=I0_mA / 1000.0,
            V_mc=V_mc,
            U0=U0,
            n_hc=n_hc,
            R_hc=R_hc,
            Q_hc=Q_hc,
        )
        optimizer = RFOptimizer(machine)
        res = optimizer.optimize_harmonic_cavity_detuning(I0_mA=I0_mA)
        df_khz = res.optimal_params["detuning_khz"]
        print(f"{df_khz:.4f},{res.xi:.4f},{res.bunch_length_factor:.2f},{res.status.name}")
    except Exception as e:
        print(f"ERROR: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
