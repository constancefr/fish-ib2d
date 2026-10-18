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
#               so editing the repo after submitting can't change a queued sweep.
#               Its code_provenance.json (+ .diff) holds the repo's git state at
#               submit time, which every run's run_info.json records
#   cases.txt   array task number -> case name
#   logs/       one Slurm log per task (<jobid>_<task>.out), plus the index job's
#   results/    one run folder per case, filled in as tasks finish; index.csv and
#               index_full.csv are built by a follow-up job once all tasks end

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
    --exclude 'results*/' --exclude viz_IB2d/ --exclude hier_IB2d_data/ \
    --exclude __pycache__/ --exclude '*.mp4' \
    cases/FinRay common external/IB2d/pyIB2d/IBM_Blackbox \
    "${SWEEP_DIR}/snapshot/")

python3 "${REPO}/common/run_records.py" provenance "${SWEEP_DIR}/snapshot" > /dev/null

(cd "${SWEEP_DIR}/snapshot/cases/FinRay" && python3 run_case.py --list) > "${SWEEP_DIR}/cases.txt"
N_CASES=$(wc -l < "${SWEEP_DIR}/cases.txt")
if [ "${N_CASES}" -eq 0 ]; then
    echo "Error: make_sweep() returned no cases." >&2
    exit 1
fi

echo "Sweep '${SWEEP_NAME}': ${N_CASES} case(s), at most ${MAX_PARALLEL} at a time"
cat "${SWEEP_DIR}/cases.txt"

ARRAY_ID=$(sbatch --parsable \
    --array="0-$((N_CASES - 1))%${MAX_PARALLEL}" \
    --job-name="finray-${SWEEP_NAME}" \
    --output="${SWEEP_DIR}/logs/%A_%a.out" \
    --export=ALL,SWEEP_DIR="${SWEEP_DIR}" \
    "${SCRIPT_DIR}/sweep_array.sbatch")
echo "Submitted array job ${ARRAY_ID}"

# Once every task has ended (finished, failed or timed out): build the sweep's
# index.csv / index_full.csv from all its run_info.json files
INDEX_ID=$(sbatch --parsable \
    --dependency="afterany:${ARRAY_ID}" \
    --partition=ICF-Free --cpus-per-task=1 --mem=1G --time=00:10:00 \
    --job-name="finray-${SWEEP_NAME}-index" \
    --output="${SWEEP_DIR}/logs/index_%j.out" \
    --wrap="source /opt/conda/etc/profile.d/conda.sh && conda activate ${CONDA_ENV} && cd ${SWEEP_DIR}/snapshot/cases/FinRay && python run_case.py index --all ${SWEEP_DIR}/results")
echo "Submitted index job ${INDEX_ID} (runs after ${ARRAY_ID})"

echo "Results will appear in ${SWEEP_DIR}/results/"
echo "Monitor with: squeue -u ${USER}"
