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

Steps 1-6 for one case are common/run_records.run_case() (shared with
NACA0020_Wing); this file only defines the sweep and the index columns.

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

EXAMPLE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = EXAMPLE_DIR / "results_approach2"

# Columns of results/index.csv, in order: the parameters and results used for
# analysis. Every other field (unused IB2d features, output switches, package
# versions, ...) stays in each run's run_info.json and in index_full.csv
# (`python run_case.py index --all`). Keep parameters here even while a sweep
# holds them constant: they are what differs between sweeps.
INDEX_COLUMNS = [
    # identity
    "folder", "run_id", "timing.start", "backfilled",
    "provenance.commit", "provenance.branch", "provenance.uncommitted_changes",
    # parameters
    "inputs.geom.ribbon_fill_density_multiplier", "inputs.geom.ribbon_fill_spacing",
    "inputs.geom.ribbon_fill_stiffness_scale", "inputs.geom.E_material",
    "inputs.geom.wall_thickness", "inputs.geom.actuation_patch_radius",
    "inputs.geom.target_stiffness_factor",
    "inputs.actuation.FREQUENCY", "inputs.actuation.AMPLITUDE",
    "inputs.actuation.STIFFNESS_RAMP_TIME",
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


def make_sweep():
    '''
    2D sweep over ribbon_fill_density_multiplier x ribbon_fill_stiffness_scale,
    everything else held fixed at BASE_CASE. Each case is named
    "mult<M>_scale<S>" and lands in its own timestamped folder
    results/mult<M>_scale<S>_<YYYYmmdd-HHMMSS>/ (see run_case()), so no run
    ever overwrites another.

    NOTE (see NOTES.md): scale=1.0 (the fill-stiffness formula's literal,
    uncalibrated value) is known to go NaN by t~3.3-3.6s at
    ribbon_fill_density_multiplier=1.0 -- higher multipliers add more,
    similarly-stiff springs and are expected to be *more* prone to this, not
    less. run_case() below flags (but does not skip) any case that went NaN,
    and records when (run_info.json, results/index.csv).
    '''
    # multipliers = [0.0, 1.0, 4.0, 8.0]
    # scales = [0.05, 0.1, 0.3, 1.0]
    # multipliers = [0.0]
    multipliers = [0.0, 1.0, 8.0, 32.0, 64.0, 128.0, 256.0]
    scales = [1.0]

    # ribbon_fill_stiffness_scale now scales EVERY spring (rails included),
    # so it matters at multiplier 0 too -- sweep it at every density.
    combos = [(mult, scale) for mult in multipliers for scale in scales]

    cases = []
    for mult, scale in combos:
        case = json.loads(json.dumps(BASE_CASE))  # deep copy of a known-clean base,
        # NOT the on-disk case_params.json -- reading that back in would also
        # pick up whatever "case_name"/one-off keys the *previous* case run
        # happened to leave behind, silently carrying them into this one.
        case["case_name"] = f"mult{mult:.2f}_scale{scale:.2f}"
        # NOTE: both must be under "geom" -- FinRay_Geom.py only reads
        # case_params.json["geom"][...]. A top-level key here is silently
        # ignored (this bit us once already: every case in a sweep would
        # quietly reuse whatever value was already baked into
        # case_params.json's "geom" section instead of the intended one).
        case["geom"]["ribbon_fill_density_multiplier"] = mult
        case["geom"]["ribbon_fill_stiffness_scale"] = scale
        cases.append(case)
    return cases


def _interrupt(signum, frame):
    raise KeyboardInterrupt


def main():
    signal.signal(signal.SIGTERM, _interrupt)  # an IDE's stop button: treat like Ctrl+C
    try:
        for case in make_sweep():
        # for case in [BASE_CASE]:    # for initial testing
            case_dir = rr.run_case(case, EXAMPLE_DIR, RESULTS_DIR, "FinRay_Geom.py")

            # Optional: chain straight into VisIt for a movie per case.
            # movie_path = case_dir / f"{case['case_name']}.mp4"
            # subprocess.run([
            #     "visit", "-cli", "-nowin", "-s", str(EXAMPLE_DIR / "make_movie.py"),
            #     "--", str(case_dir), str(movie_path),
            # ], check=True)
    finally:
        # rebuilt even if a case crashed, so finished runs are always listed
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