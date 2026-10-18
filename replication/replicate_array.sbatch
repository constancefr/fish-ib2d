#!/bin/bash
# SLURM job array: one task per run in manifest.json (task i = run i of `replicate.py list`).
# Submit from the repo root:
#     mkdir -p logs
#     sbatch replication/replicate_array.sbatch                      # all 32 runs
#     sbatch --array=0-5 replication/replicate_array.sbatch          # e.g. only approach1 (rows 0-5)
# Finished runs are skipped, so resubmitting only redoes runs that are missing or were
# stopped early (e.g. by the time limit). EXTRA_ARGS=--force reruns finished ones too.
#
# Adapt the #SBATCH lines and the environment setup below to your cluster.
#SBATCH --job-name=fish-replicate
#SBATCH --array=0-31%8
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%A_%a.out

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"

# --- environment (edit) ---
# module load miniconda
# source activate ib2d
export MPLBACKEND=Agg OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

OUT="${OUT:-results_replicated}"   # e.g. OUT=/scratch/$USER/fish-ib2d sbatch ...
python replication/replicate.py run --all --index "${SLURM_ARRAY_TASK_ID}" --out "${OUT}" ${EXTRA_ARGS:-}
