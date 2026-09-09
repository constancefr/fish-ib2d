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

import itertools
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

EXAMPLE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = EXAMPLE_DIR / "results"


# --- Define cases here
# Every key under "geom"/"actuation" becomes available to FinRay_Geom.py /
#   update_Target_Point_Positions.py via case_params.json.
# Keys under "input2d" are patched directly into the input2d text file.
BASE_CASE = {
    "case_name": "baseline",
    "input2d": {
        "mu": 0.001,
        "rho": 1.0,
        "Tfinal": 5.0,
        "dt": 1.0e-3,
        "Nx": 32, # keep even
        "Ny": 32, # keep even
        "Lx": 1.0,
        "Ly": 1.0,
        "supp": 4, # keep even
    },
    "geom": {
        "E_material": 0.74e6,
        "wall_thickness": 0.00396,
        "extrude_depth": 0.05,
    },
    "actuation": {
        "FREQUENCY": 1.0,
        "AMPLITUDE": 0.01,
        "TARGET_STIFFNESS_SCALE": 0.2,
        "STIFFNESS_RAMP_TIME": 0.2,
    },
}


def make_sweep():
    '''
    Example sweep over AMPLITUDE x E_material, everything else held fixed.
    '''
    cases = []
    for amp, E in itertools.product([0.005, 0.01, 0.02], [0.5e6, 0.74e6, 1.0e6]):
        case = json.loads(json.dumps(BASE_CASE))  # deep copy
        case["case_name"] = f"amp{amp:.3f}_E{E:.2e}"
        case["actuation"]["AMPLITUDE"] = amp
        case["geom"]["E_material"] = E
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

    print(f"=== Done: results in {case_dir} ===")
    return case_dir


def main():
    # for case in make_sweep():
    for case in [BASE_CASE]:    # for initial testing
        case_dir = run_case(case)

        # Optional: chain straight into VisIt for a movie per case.
        # movie_path = case_dir / f"{case['case_name']}.mp4"
        # subprocess.run([
        #     "visit", "-cli", "-nowin", "-s", str(EXAMPLE_DIR / "make_movie.py"),
        #     "--", str(case_dir), str(movie_path),
        # ], check=True)


if __name__ == "__main__":
    main()