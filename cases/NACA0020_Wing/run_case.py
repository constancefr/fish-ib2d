'''
Parameter-sweep driver for the NACA0020 wing (same workflow as
FinRay/run_case.py). For each case, common/run_records.run_case():

  1. writes case_params.json and patches the case's "input2d" values into
     input2d
  2. runs Wing_Geom.py, then main2d.py, each in a fresh subprocess
  3. moves viz_IB2d/ and hier_IB2d_data/ into a new folder
     results/<case_name>_<YYYYmmdd-HHMMSS>/, with copies of case_params.json,
     input2d, wing.vertex and wing.target (analyze_run.py reads the last two)
  4. writes that folder's run_info.json (inputs, mesh metrics, timing,
     outcome, code version), then runs analyze_run.py on it

After the sweep, results/index.csv is rebuilt (one row per run).

NOTE: stiffness_scale is calibrated for one mesh (verify_structure.py prints
it). A sweep that changes fill_density_multiplier, fill_spacing or the grid
must set the matching stiffness_scale for each case.

Run from this folder:  python run_case.py
                       python run_case.py index [--all] [DIR ...]   (rebuild the index only)
'''

import json
import signal
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
import run_records as rr

EXAMPLE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = EXAMPLE_DIR / "results"

# Columns of results/index.csv, in order (every field is in run_info.json and
# in index_full.csv: `python run_case.py index --all`).
INDEX_COLUMNS = [
    # identity
    "folder", "run_id", "timing.start", "backfilled",
    "provenance.commit", "provenance.branch", "provenance.uncommitted_changes",
    # parameters
    "inputs.geom.fill_density_multiplier", "inputs.geom.fill_spacing",
    "inputs.geom.stiffness_scale", "inputs.geom.E_flex", "inputs.geom.rigid_length",
    "inputs.geom.init_q", "inputs.geom.target_stiffness_factor",
    "inputs.flow.U", "inputs.flow.aoa_deg", "inputs.flow.aoa_switch_time",
    "inputs.flow.fringe_rate", "inputs.flow.mean_rate",
    "inputs.input2d.mu", "inputs.input2d.rho", "inputs.input2d.dt", "inputs.input2d.Tfinal",
    "inputs.input2d.Nx", "inputs.input2d.Ny", "inputs.input2d.Lx", "inputs.input2d.Ly",
    "inputs.input2d.supp",
    # mesh
    "mesh.edge_spacing", "mesh.n_edge_points", "mesh.n_fill_points", "mesh.n_points",
    "mesh.n_springs", "mesh.k_median", "mesh.k_max", "mesh.n_target_points", "mesh.k_Target",
    "mesh.mesh_area_fraction", "mesh.build_time_s",
    # timing
    "timing.geometry_wall_time_s", "timing.sim_wall_time_s", "timing.sim_time_per_step_ms",
    # outcome
    "outcome.status", "outcome.first_nan_time", "outcome.last_output_time",
    "outcome.interrupted",
]

# Base case = case_params.json on disk: editing that file IS how the sweep's
# starting point changes.
BASE_CASE = json.loads((EXAMPLE_DIR / "case_params.json").read_text())


def make_sweep():
    '''Example: viscosity multiple (Re). Edit freely; every key under
    "geom"/"flow" reaches Wing_Geom.py / the forcing via case_params.json.'''
    cases = []
    for mu_scale in (30, 100):
        case = json.loads(json.dumps(BASE_CASE))      # deep copy
        case["case_name"] = f"mu{mu_scale}x"
        case["input2d"]["mu"] = 1.0023e-3 * mu_scale
        cases.append(case)
    return cases


def _interrupt(signum, frame):
    raise KeyboardInterrupt


def main():
    signal.signal(signal.SIGTERM, _interrupt)  # an IDE's stop button: treat like Ctrl+C
    try:
        for case in [BASE_CASE]:      # swap for make_sweep() to run the sweep
            rr.run_case(case, EXAMPLE_DIR, RESULTS_DIR, "Wing_Geom.py",
                        post_scripts=("analyze_run.py",), copy_files=("wing.vertex", "wing.target"))
    finally:
        # rebuilt even if a case crashed, so finished runs are always listed
        if RESULTS_DIR.exists():
            print(f"\nIndex: {rr.write_index(RESULTS_DIR, INDEX_COLUMNS)}")


if __name__ == "__main__":
    if sys.argv[1:2] == ["index"]:
        dirs = [a for a in sys.argv[2:] if a != "--all"] or [RESULTS_DIR]
        for d in dirs:
            print("wrote", rr.write_index(d, INDEX_COLUMNS))
            if "--all" in sys.argv:
                print("wrote", rr.write_index(d))
    else:
        main()
