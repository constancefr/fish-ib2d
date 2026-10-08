'''
Parameter-sweep driver for the NACA0020 wing example (same workflow as
FinRay/run_case.py). For each case: write case_params.json, patch input2d,
regenerate the structure (Wing_Geom.py), run main2d.py in a fresh subprocess,
run analyze_run.py, and archive viz_IB2d / hier_IB2d_data / history.* under
results/<case_name>/.

Run from this folder:  python run_case.py
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

BASE_CASE = json.loads((EXAMPLE_DIR / "case_params.json").read_text())


def make_sweep():
    '''Example: viscosity multiple (Re) x flexural modulus. Edit freely.'''
    cases = []
    for mu_scale in (30, 100):
        for E in (1.0e4, 3.0e4):
            case = json.loads(json.dumps(BASE_CASE))
            case["case_name"] = f"mu{mu_scale}x_E{E:.0e}"
            case["input2d"]["mu"] = 1.0023e-3 * mu_scale
            case["geom"]["E_flex"] = E
            cases.append(case)
    return cases


def patch_input2d(values: dict):
    path = EXAMPLE_DIR / "input2d"
    text = path.read_text()
    for key, val in values.items():
        pattern = re.compile(rf"^(\s*{re.escape(key)}\s*=\s*)([^\s%]+)", re.MULTILINE)
        if not pattern.search(text):
            raise ValueError(f"Couldn't find '{key}' in input2d")
        text = pattern.sub(lambda m: m.group(1) + repr(val), text)
    path.write_text(text)


def run_case(case: dict):
    print(f"\n=== Running case: {case['case_name']} ===")
    params_path = EXAMPLE_DIR / "case_params.json"
    params_path.write_text(json.dumps(case, indent=2))
    patch_input2d(case["input2d"])

    env = dict(os.environ, MPLBACKEND="Agg")
    subprocess.run([sys.executable, "Wing_Geom.py"], cwd=EXAMPLE_DIR, env=env, check=True)
    subprocess.run([sys.executable, "main2d.py"], cwd=EXAMPLE_DIR, check=True)
    subprocess.run([sys.executable, "analyze_run.py"], cwd=EXAMPLE_DIR, env=env, check=True)

    case_dir = RESULTS_DIR / case["case_name"]
    case_dir.mkdir(parents=True, exist_ok=True)
    for name in ("viz_IB2d", "hier_IB2d_data"):
        src = EXAMPLE_DIR / name
        if src.exists():
            dest = case_dir / name
            if dest.exists():
                shutil.rmtree(dest)
            shutil.move(str(src), str(dest))
    for name in ("history.npz", "history.png", "input2d", "case_params.json"):
        if (EXAMPLE_DIR / name).exists():
            shutil.copy(EXAMPLE_DIR / name, case_dir / name)
    print(f"=== Done: results in {case_dir} ===")
    return case_dir


def main():
    for case in [BASE_CASE]:      # swap for make_sweep() to run the sweep
        run_case(case)


if __name__ == "__main__":
    main()
