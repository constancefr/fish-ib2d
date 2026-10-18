#!/bin/bash
# Submit run_case.py's make_sweep() as a Slurm job array, one task per case.
#
# Usage (on the head node, from anywhere):
#   cases/FinRay/slurm/submit_sweep.sh <sweep_name>
#
# Optional environment variables:
#   MAX_PARALLEL  cases allowed to run at once (default 8). Change a submitted
#                 array with: scontrol update JobId=<jobid> ArrayTaskThrottle=<N>
#   CONDA_ENV     conda env the jobs activate (default ib2d)
#   DATA_ROOT     where sweeps are stored (default ~/sim_data/fish-ib2d/FinRay)
#
# Creates $DATA_ROOT/<sweep_name>/ containing:
#   snapshot/   frozen copy of the code at submit time; every task runs from this,
#               so editing the repo after submitting can't change a queued sweep
#   git.txt     commit, status and diff of the repo at submit time
#   cases.txt   array task index -> case name
#   logs/       one Slurm log per task (<jobid>_<task>.out)
#   results/    one folder per case, filled in as tasks finish

set -euo pipefail

if [ $# -ne 1 ]; then
    sed -n '2,22p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
fi
SWEEP_NAME="$1"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
DATA_ROOT="${DATA_ROOT:-${HOME}/sim_data/fish-ib2d/FinRay}"
SWEEP_DIR="${DATA_ROOT}/${SWEEP_NAME}"
MAX_PARALLEL="${MAX_PARALLEL:-8}"
export CONDA_ENV="${CONDA_ENV:-ib2d}"

if [ -e "${SWEEP_DIR}" ]; then
    echo "Error: ${SWEEP_DIR} already exists -- pick another sweep name." >&2
    exit 1
fi
if [ ! -d "${REPO}/external/IB2d/pyIB2d/IBM_Blackbox" ]; then
    echo "Error: IB2d submodule missing -- run: git -C ${REPO} submodule update --init" >&2
    exit 1
fi

# Snapshot only what a FinRay run needs (not .git, other cases, or old output)
mkdir -p "${SWEEP_DIR}/snapshot" "${SWEEP_DIR}/logs" "${SWEEP_DIR}/results"
(cd "${REPO}" && rsync -aR \
    --exclude results/ --exclude viz_IB2d/ --exclude hier_IB2d_data/ \
    --exclude __pycache__/ --exclude '*.mp4' \
    cases/FinRay common external/IB2d/pyIB2d/IBM_Blackbox \
    "${SWEEP_DIR}/snapshot/")

{
    echo "commit: $(git -C "${REPO}" rev-parse HEAD)"
    echo; echo "status:"; git -C "${REPO}" status --short
    echo; echo "diff:";   git -C "${REPO}" diff
} > "${SWEEP_DIR}/git.txt"

(cd "${SWEEP_DIR}/snapshot/cases/FinRay" && python3 run_case.py --list) > "${SWEEP_DIR}/cases.txt"
N_CASES=$(wc -l < "${SWEEP_DIR}/cases.txt")
if [ "${N_CASES}" -eq 0 ]; then
    echo "Error: make_sweep() returned no cases." >&2
    exit 1
fi

echo "Sweep '${SWEEP_NAME}': ${N_CASES} case(s), at most ${MAX_PARALLEL} at a time"
cat "${SWEEP_DIR}/cases.txt"

sbatch \
    --array="0-$((N_CASES - 1))%${MAX_PARALLEL}" \
    --job-name="finray-${SWEEP_NAME}" \
    --output="${SWEEP_DIR}/logs/%A_%a.out" \
    --export=ALL,SWEEP_DIR="${SWEEP_DIR}" \
    "${SCRIPT_DIR}/sweep_array.sbatch"

echo "Results will appear in ${SWEEP_DIR}/results/"
echo "Monitor with: squeue -u ${USER}"
