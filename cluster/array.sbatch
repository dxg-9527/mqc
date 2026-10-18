#!/usr/bin/env bash
set -euo pipefail
# All resource options are supplied by submit.py, not hard-coded per model.
MQC_PYTHON="$1"
MQC_REPO="$2"
MQC_PLAN="$3"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONHASHSEED=0
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MPLBACKEND=Agg
cd "$MQC_REPO"
exec srun --cpu-bind=cores "$MQC_PYTHON" cluster/worker.py \
    --plan "$MQC_PLAN" --slot "${SLURM_ARRAY_TASK_ID:?Missing Slurm array index}"
