#!/bin/bash
#MSUB -r STABLE_MultiA100
#MSUB -q a100
#MSUB -A gen1234
#MSUB -n 1
#MSUB -c 128
#MSUB -T 7200
#MSUB -m scratch,work
#MSUB -o stable_a100_%I.out
#MSUB -e stable_a100_%I.err

# ==============================================================================
# TGCC Irene A100 Multi-GPU Submission Template for STABLE_Py_TGCC
# Submit via:
#   ccc_msub scripts/submit_tgcc_a100.sh
# ==============================================================================

set -e
module purge
module load python cuda

cd "${CCCSCRATCHDIR}/STABLE_Py_TGCC"
echo "[TGCC Irene A100] Node: $(hostname), Date: $(date)"
nvidia-smi

# Start the RPC streaming server on the allocated A100 node (or run a batch scan request)
if [ -f "job_request.json" ]; then
    python3 -m stable_tgcc.server.ccc_scheduler \
        --job-file "job_request.json" \
        --telemetry-file "telemetry.jsonl"
else
    python3 run_server.py --host 0.0.0.0 --port 18765
fi
