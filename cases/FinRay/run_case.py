'''
Parameter-sweep driver for the FinRay IB2d example. Configure one or more
"cases" (fluid, material, and actuation parameters), and this script will:

  1. write a case_params.json that FinRay_Geom.py and
     update_Target_Point_Positions.py read at the top of their files
  2. patch the relevant lines in input2d (fluid viscosity/density,
     Tfinal, dt, ...) via regex text substitution -- comments and
     ordering in input2d are left untouched
  3. run `python FinRay_Geom.py` in a fresh subprocess to regenerate
     .vertex/.spring/.target/.geo_connect
  4. run `python main2d.py` in a fresh subprocess to run the simulation
  5. move the resulting viz_IB2d/ (and hier_IB2d_data/, if present) into
     a new folder results/<case_name>_<YYYYmmdd-HHMMSS>/, alongside a copy
     of case_params.json and input2d
  6. write that folder's run_info.json (inputs, mesh metrics, timing,
     outcome, code version), and after the sweep rebuild results/index.csv,
     one row per run
  7. restore case_params.json and input2d to what they were before the
     case, so a run leaves the tracked inputs unchanged

Steps 1-7 for one case are common/run_records.run_case() (shared with
NACA0020_Wing). After each run, postprocess.py adds the tip and fluid
measures to its run_info.json. The sweeps themselves are defined in
sweeps.py and chosen with --sweep; this file runs them and sets the index
columns.

To run one case per Slurm array task on the cluster, see slurm/submit_sweep.sh
(`--list` and `--case N` below).

WHY SUBPROCESSES ---
update_Target_Point_Positions.py caches its baseline corner positions in
module-level globals the first time it's called (`if _baseline['x'] is
None: ...`). If main2d() were called repeatedly inside one long-lived
Python process, every case after the first would silently reuse case
#1's baseline geometry instead of its own. Running each case as its own
subprocess gives it a clean interpreter, mirroring what would happen when 
running `python main2d.py` by hand for each case.

Run this script from inside the wanted Example folder, e.g. FinRay/
'''

import json
import signal
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
import run_records as rr
import sweeps

EXAMPLE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = EXAMPLE_DIR / "results_approach2"

# Columns of results/index.csv, in order: the parameters and results used for
# analysis. Every other field (unused IB2d features, output switches, package
# versions, ...) stays in each run's run_info.json and in index_full.csv
# (`python run_case.py index --all`). Keep parameters here even while a sweep
# holds them constant: they are what differs between sweeps.
INDEX_COLUMNS = [
    # identity
    "folder", "run_id", "inputs.sweep", "timing.start", "backfilled",
    "provenance.commit", "provenance.branch", "provenance.uncommitted_changes",
    # parameters
    "inputs.geom.ribbon_fill_density_multiplier", "inputs.geom.ribbon_fill_spacing",
    "inputs.geom.halton_seed",
    "inputs.geom.ribbon_fill_stiffness_scale", "inputs.geom.E_material",
    "inputs.geom.wall_thickness", "inputs.geom.actuation_patch_radius",
    "inputs.geom.target_stiffness_factor",
    "inputs.actuation.FREQUENCY", "inputs.actuation.AMPLITUDE",
    "inputs.actuation.STIFFNESS_RAMP_TIME",
    "inputs.run_length.transient_cycles", "inputs.run_length.measured_cycles",
    "inputs.input2d.mu", "inputs.input2d.rho", "inputs.input2d.dt", "inputs.input2d.Tfinal",
    "inputs.input2d.Nx", "inputs.input2d.Ny", "inputs.input2d.Lx", "inputs.input2d.Ly",
    "inputs.input2d.supp",
    # mesh
    "mesh.edge_spacing", "mesh.n_edge_points", "mesh.n_fill_points_requested",
    "mesh.n_fill_points", "mesh.n_points", "mesh.n_springs", "mesh.k_median", "mesh.k_max",
    "mesh.n_target_points", "mesh.k_Target", "mesh.mesh_area_fraction", "mesh.build_time_s",
    # timing (backfilled runs only have the approximate wall time)
    "timing.geometry_wall_time_s", "timing.sim_wall_time_s", "timing.sim_wall_time_s_approx",
    "timing.sim_time_per_step_ms",
    # outcome
    "outcome.status", "outcome.first_nan_time", "outcome.last_output_time",
    "outcome.interrupted",
    # analysis (postprocess.py)
    "analysis.tip.n_cycles", "analysis.tip.y_amplitude", "analysis.tip.y_amplitude_cv",
    "analysis.tip.y_mean", "analysis.tip.x_amplitude", "analysis.tip.phase_lag_deg",
    "analysis.fluid.gs_ke_growth_rate", "analysis.fluid.gs_ke_fraction_end",
    "analysis.fluid.u_max", "analysis.fluid.cell_re_max",
]


# --- Base case: loaded from case_params.json on disk ---
# Every key under "geom"/"actuation" becomes available to FinRay_Geom.py /
#   update_Target_Point_Positions.py via case_params.json.
# Keys under "input2d" are patched directly into the input2d text file.
#
# BASE_CASE used to be a hardcoded dict here, separate from case_params.json
# -- which meant hand-editing case_params.json (e.g. to try a value, or
# because that's what a previous run_case() call last wrote) had no effect
# on the next sweep: it would silently start over from whatever was baked
# into this file instead. Loading it from case_params.json means editing
# that file directly IS how you change the sweep's starting point; a
# `git diff case_params.json` also then shows exactly what changed.
#
# NOTE (see NOTES.md): ribbon_fill_stiffness_scale=1.0 (the fill-stiffness
# formula's literal, uncalibrated value) is known to go NaN by t~3.3-3.6s at
# ribbon_fill_density_multiplier=1.0, dt=1e-3, E_material=0.74e6; 0.1 stayed
# clean through a full Tfinal=5.0 run at those settings. Re-check (rerun to
# Tfinal and grep the last viz_IB2d/lagPtsConnect.*.vtk frame for "nan", or
# just watch for run_case()'s own divergence warning below) whenever
# case_params.json's E_material, dt, or ribbon_fill_density_multiplier
# change materially from those values.
def _load_base_case():
    with open(EXAMPLE_DIR / "case_params.json") as f:
        return json.load(f)


BASE_CASE = _load_base_case()


POST_SCRIPTS = ("postprocess.py",)


def _interrupt(signum, frame):
    raise KeyboardInterrupt


def main():
    '''
    No arguments: run every case of the sweep one after another, then
                  rebuild the index (local use).
    --sweep NAME: the sweep in sweeps.py to run (default: base, the base case
                  alone)
    --list:       print "<N> <case_name>" for each case and exit, without
                  running anything -- slurm/submit_sweep.sh uses this to size
                  the job array.
    --case N:     run only case N of the sweep (0-based), without rebuilding
                  the index -- one Slurm array task per case, see
                  slurm/sweep_array.sbatch; submit_sweep.sh queues a follow-up
                  job that builds the sweep's index once every case is done.
    '''
    import argparse
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--sweep", default="base", metavar="NAME",
                        help=f"sweep to run, from sweeps.py: {', '.join(sweeps.SWEEPS)}")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--list", action="store_true",
                       help="print the number and name of each case, then exit")
    group.add_argument("--case", type=int, metavar="N",
                       help="run only case N of the sweep (0-based)")
    args = parser.parse_args()

    cases = sweeps.build(args.sweep, BASE_CASE)

    if args.list:
        for i, case in enumerate(cases):
            print(i, case["case_name"])
        return
    if args.case is not None:
        if not 0 <= args.case < len(cases):
            sys.exit(f"--case {args.case} out of range: the sweep has "
                     f"{len(cases)} case(s), numbered 0-{len(cases) - 1}")
        cases = [cases[args.case]]

    signal.signal(signal.SIGTERM, _interrupt)  # an IDE's stop button (or Slurm): treat like Ctrl+C
    try:
        for case in cases:
            case_dir = rr.run_case(case, EXAMPLE_DIR, RESULTS_DIR, "FinRay_Geom.py",
                                   post_scripts=POST_SCRIPTS)

            # Optional: chain straight into VisIt for a movie per case.
            # movie_path = case_dir / f"{case['case_name']}.mp4"
            # subprocess.run([
            #     "visit", "-cli", "-nowin", "-s", str(EXAMPLE_DIR / "make_movie.py"),
            #     "--", str(case_dir), str(movie_path),
            # ], check=True)
    finally:
        # rebuilt even if a case crashed, so finished runs are always listed
        # (not per array task: each would hold only its own run)
        if args.case is None:
            print(f"\nIndex: {rr.write_index(RESULTS_DIR, INDEX_COLUMNS)}")


if __name__ == "__main__":
    # `python run_case.py index [--all] [DIR ...]` only rebuilds index.csv
    # (and index_full.csv with --all) for DIR (default RESULTS_DIR)
    if sys.argv[1:2] == ["index"]:
        dirs = [a for a in sys.argv[2:] if a != "--all"] or [RESULTS_DIR]
        for d in dirs:
            print("wrote", rr.write_index(d, INDEX_COLUMNS))
            if "--all" in sys.argv:
                print("wrote", rr.write_index(d))
    else:
        main()