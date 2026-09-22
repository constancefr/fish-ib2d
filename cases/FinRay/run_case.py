'''
Parameter-sweep driver for the FinRay IB2d example. Configure one or more
"cases" (fluid, material, and actuation parameters), and this script will:

  1. write a case_params.json that FinRay_Geom.py and
     update_Target_Point_Positions.py read at the top of their files
  2. patch the relevant lines in input2d (fluid viscosity/density,
     Tfinal, dt, ...) via regex text substitution -- comments and
     ordering in input2d are left untouched
  3. run `python FinRay_Geom.py` in a fresh subprocess to regenerate
     .vertex/.spring/.beam/.target/.geo_connect
  4. run `python main2d.py` in a fresh subprocess to run the simulation
  5. move the resulting viz_IB2d/ (and hier_IB2d_data/, if present) into
     a per-case results folder, alongside a copy of case_params.json and
     input2d, for provenance

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
from pathlib import Path

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
    "mult<M>_scale<S>" and lands in results/mult<M>_scale<S>/ (see run_case()),
    so results from different (multiplier, scale) combinations never collide
    or overwrite each other.

    NOTE (see NOTES.md): scale=1.0 (the fill-stiffness formula's literal,
    uncalibrated value) is known to go NaN by t~3.3-3.6s at
    ribbon_fill_density_multiplier=1.0 -- higher multipliers add more,
    similarly-stiff springs and are expected to be *more* prone to this, not
    less. run_case() below flags (but does not skip) any case whose final
    output frame contains NaN, so a blown-up case is obvious without having
    to manually inspect each results/ folder.
    '''
    # multipliers = [0.0, 1.0, 4.0, 8.0]
    # scales = [0.05, 0.1, 0.3, 1.0]
    multipliers = [0.0, 1.0, 8.0, 32.0, 64.0]
    scales = [1.0]

    combos = []
    for mult in multipliers:
        if mult == 0.0:
            # ribbon_fill_stiffness_scale is meaningless with no fill points
            # (no 'fill' springs get created at all -- see
            # build_Tail_Ribbon_Connections) -- sweeping scale here would
            # just rerun the identical no-fill baseline 4x.
            combos.append((mult, scales[0]))
        else:
            combos.extend((mult, scale) for scale in scales)

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
    print(f"\n=== Running case: {case['case_name']} ===")

    # a) write case_params.json (read by FinRay_Geom.py and
    #    update_Target_Point_Positions.py)
    params_path = EXAMPLE_DIR / "case_params.json"
    params_path.write_text(json.dumps(case, indent=2))

    # b) patch input2d
    patch_input2d(case["input2d"])

    # c) regenerate geometry/input files
    env = dict(os.environ, MPLBACKEND="Agg") # Agg backend so that plt.show() is a no-op and doesn't block
    subprocess.run([sys.executable, "FinRay_Geom.py"], cwd=EXAMPLE_DIR, env=env, check=True)

    # d) run the simulation (fresh subprocess)
    subprocess.run([sys.executable, "main2d.py"], cwd=EXAMPLE_DIR, check=True)

    # e) archive results
    case_dir = RESULTS_DIR / case["case_name"]
    case_dir.mkdir(parents=True, exist_ok=True)
    for folder in ("viz_IB2d", "hier_IB2d_data"):
        src = EXAMPLE_DIR / folder
        if src.exists():
            dest = case_dir / folder
            if dest.exists():
                shutil.rmtree(dest)
            shutil.move(str(src), str(dest))
    shutil.copy(params_path, case_dir / "case_params.json")
    shutil.copy(EXAMPLE_DIR / "input2d", case_dir / "input2d")

    if _final_frame_has_nan(case_dir):
        print(f"=== WARNING: {case['case_name']} diverged (NaN in its final "
              f"output frame) -- results in {case_dir} are not usable ===")
    else:
        print(f"=== Done: results in {case_dir} ===")
    return case_dir


def _final_frame_has_nan(case_dir):
    '''
    Cheap divergence check: does the *last* lagPtsConnect.*.vtk frame
    contain "nan"? A spring-force divide-by-zero (two points colliding)
    poisons position/velocity data from that point on but does not raise an
    exception or a nonzero exit code, so a diverged run otherwise looks
    identical to a healthy one until you open the output by hand -- see
    NOTES.md for how this was first found (ribbon_fill_stiffness_scale=1.0
    going NaN partway through a run).
    '''
    frames = sorted((case_dir / "viz_IB2d").glob("lagPtsConnect.*.vtk"))
    if not frames:
        return False
    return "nan" in frames[-1].read_text().lower()


def main():
    for case in make_sweep():
    # for case in [BASE_CASE]:    # for initial testing
        case_dir = run_case(case)

        # Optional: chain straight into VisIt for a movie per case.
        # movie_path = case_dir / f"{case['case_name']}.mp4"
        # subprocess.run([
        #     "visit", "-cli", "-nowin", "-s", str(EXAMPLE_DIR / "make_movie.py"),
        #     "--", str(case_dir), str(movie_path),
        # ], check=True)


if __name__ == "__main__":
    main()