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
     outcome, code version; see common/run_records.py), and after the
     sweep rebuild results/index.csv, one row per run

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
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
import run_records as rr

EXAMPLE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = EXAMPLE_DIR / "results"


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
    multipliers = [1.0]
    # multipliers = [0.0, 1.0, 8.0, 32.0, 64.0]
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


# --- input2d patching
def patch_input2d(values: dict):
    '''
    Rewrites matching `name = value` lines in input2d in place.
    '''
    path = EXAMPLE_DIR / "input2d"
    text = path.read_text()

    for key, val in values.items():
        # Matches "key = 0.01" or "key= 0.01   % comment"; keeps the comment.
        pattern = re.compile(rf"^(\s*{re.escape(key)}\s*=\s*)([^\s%]+)", re.MULTILINE)
        if not pattern.search(text):
            raise ValueError(
                f"Couldn't find '{key}' in input2d -- check the exact "
                f"parameter name/formatting in your input2d file."
            )
        text = pattern.sub(lambda m: m.group(1) + repr(val), text)

    path.write_text(text)


# --- Run one case
def run_case(case: dict):
    '''
    Runs one case into its own new folder results/<case_name>_<YYYYmmdd-HHMMSS>/
    (never overwrites an earlier run) and writes its run_info.json there:
    inputs, mesh metrics, timing, outcome, and the code version (commit,
    plus code.diff if the tree had uncommitted changes). See
    common/run_records.py.
    '''
    start = datetime.now()
    run_id = f"{case['case_name']}_{start:%Y%m%d-%H%M%S}"
    case_dir = RESULTS_DIR / run_id
    case_dir.mkdir(parents=True)
    print(f"\n=== Running case: {run_id} ===")
    # code version as it is at launch (edits made while a sweep runs don't count)
    code = rr.provenance(EXAMPLE_DIR, case_dir / "code.diff")

    # a) write case_params.json (read by FinRay_Geom.py and
    #    update_Target_Point_Positions.py)
    params_path = EXAMPLE_DIR / "case_params.json"
    params_path.write_text(json.dumps(case, indent=2))

    # b) patch input2d
    patch_input2d(case["input2d"])

    # c) regenerate geometry/input files (also writes mesh_info.json)
    env = dict(os.environ, MPLBACKEND="Agg") # Agg backend so that plt.show() is a no-op and doesn't block
    t0 = time.perf_counter()
    subprocess.run([sys.executable, "FinRay_Geom.py"], cwd=EXAMPLE_DIR, env=env, check=True)
    t_geom = time.perf_counter() - t0

    # d) run the simulation (fresh subprocess)
    t0 = time.perf_counter()
    sim = subprocess.run([sys.executable, "main2d.py"], cwd=EXAMPLE_DIR)
    t_sim = time.perf_counter() - t0

    # e) archive results
    for folder in ("viz_IB2d", "hier_IB2d_data"):
        src = EXAMPLE_DIR / folder
        if src.exists():
            shutil.move(str(src), str(case_dir / folder))
    shutil.copy(params_path, case_dir / "case_params.json")
    shutil.copy(EXAMPLE_DIR / "input2d", case_dir / "input2d")

    # f) run record
    inp = rr.read_input2d(case_dir / "input2d")
    n_steps = round(inp["Tfinal"] / inp["dt"])
    info = dict(
        run_id=run_id,
        case_name=case["case_name"],
        inputs=rr.inputs(case, inp),
        mesh=json.loads((EXAMPLE_DIR / "mesh_info.json").read_text()),
        timing=dict(start=start.isoformat(timespec="seconds"),
                    end=datetime.now().isoformat(timespec="seconds"),
                    geometry_wall_time_s=round(t_geom, 1),
                    sim_wall_time_s=round(t_sim, 1),
                    n_steps=n_steps,
                    sim_time_per_step_ms=round(1e3 * t_sim / n_steps, 2)),
        outcome=dict(rr.outcome(case_dir / "viz_IB2d", inp["dt"], inp["print_dump"], inp["Tfinal"]),
                     exit_code=sim.returncode),
        provenance=code,
    )
    rr.write_run_info(case_dir, info)

    status = info["outcome"]["status"]
    if status == "diverged":
        print(f"=== WARNING: {run_id} diverged (NaN from t = "
              f"{info['outcome']['first_nan_time']:.2f} s) -- see {case_dir} ===")
    else:
        print(f"=== {status}: results in {case_dir} ===")
    if sim.returncode:
        raise subprocess.CalledProcessError(sim.returncode, "main2d.py")
    return case_dir


def main():
    try:
        for case in make_sweep():
        # for case in [BASE_CASE]:    # for initial testing
            case_dir = run_case(case)

            # Optional: chain straight into VisIt for a movie per case.
            # movie_path = case_dir / f"{case['case_name']}.mp4"
            # subprocess.run([
            #     "visit", "-cli", "-nowin", "-s", str(EXAMPLE_DIR / "make_movie.py"),
            #     "--", str(case_dir), str(movie_path),
            # ], check=True)
    finally:
        # rebuilt even if a case crashed, so finished runs are always listed
        print(f"\nIndex: {rr.write_index(RESULTS_DIR)}")


if __name__ == "__main__":
    main()