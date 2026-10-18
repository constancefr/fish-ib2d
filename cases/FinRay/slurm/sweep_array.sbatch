#!/bin/bash
# One Slurm array task = one case of run_case.py's make_sweep().
# Don't sbatch this directly -- use submit_sweep.sh, which snapshots the code,
# sizes the array and sets SWEEP_DIR.
#
# Each task:
#   1. copies the sweep's code snapshot to the compute node's local /disk/scratch
#   2. runs `python run_case.py --case $SLURM_ARRAY_TASK_ID` there, so all the
#      .vtk output is written to local disk, not the shared filesystem
#   3. copies the case's run folder back to $SWEEP_DIR/results/ and deletes the
#      scratch copy -- also on failure or timeout (run_case.py records a
#      stopped run as interrupted), so partial output and scratch cleanup are
#      never lost
#
#SBATCH --partition=ICF-Free
#SBATCH --cpus-per-task=1          # IB2d is single-threaded
#SBATCH --mem=4G                   # check real usage: sacct -j <jobid> --format=JobID,Elapsed,MaxRSS
#SBATCH --time=02:00:00
#SBATCH --signal=B:TERM@300        # warn the script 5 min before the time limit, to copy back
# --array, --job-name and --output are set by submit_sweep.sh

set -euo pipefail

: "${SWEEP_DIR:?SWEEP_DIR not set -- submit with submit_sweep.sh}"
CONDA_ENV="${CONDA_ENV:-ib2d}"
TASK="${SLURM_ARRAY_TASK_ID}"
WORK="/disk/scratch/${USER}/finray_${SLURM_ARRAY_JOB_ID}_${TASK}"
CASE_DIR="${WORK}/cases/FinRay"

echo "Task ${TASK} on $(hostname) at $(date), scratch: ${WORK}"

copy_back() {
    local rc=$?
    trap - EXIT TERM
    set +e
    # On timeout the simulation is still running: stop it before copying its files
    if [ -n "${PY_PID:-}" ] && kill -0 "${PY_PID}" 2>/dev/null; then
        kill "${PY_PID}"
        wait "${PY_PID}"
    fi
    # run_case.py's RESULTS_DIR (results*/ isn't in the snapshot, so this is only this run)
    for d in "${CASE_DIR}"/results*/; do
        [ -d "${d}" ] && rsync -a "${d}" "${SWEEP_DIR}/results/"
    done
    # run_case.py was killed outright, so never archived its output: keep
    # whatever frames it wrote for inspection
    if [ -d "${CASE_DIR}/viz_IB2d" ]; then
        local dest="${SWEEP_DIR}/results/_incomplete_task${TASK}"
        mkdir -p "${dest}"
        rsync -a "${CASE_DIR}/viz_IB2d" "${CASE_DIR}/case_params.json" "${CASE_DIR}/input2d" "${dest}/"
        echo "Run did not finish: partial output copied to ${dest}"
    fi
    rm -rf "${WORK}"
    echo "Task ${TASK} finished at $(date) with exit code ${rc}"
    exit "${rc}"
}
trap copy_back EXIT
trap 'exit 143' TERM

mkdir -p "${WORK}"
rsync -a "${SWEEP_DIR}/snapshot/" "${WORK}/"

source /opt/conda/etc/profile.d/conda.sh
conda activate "${CONDA_ENV}"

# Keep numpy/numba to the one core we asked for
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export MPLBACKEND=Agg

cd "${CASE_DIR}"
# Run in the background and wait, so the TERM trap fires straight away on timeout
python run_case.py --case "${TASK}" &
PY_PID=$!
wait "${PY_PID}"
