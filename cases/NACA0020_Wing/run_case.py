'''
Parameter-sweep driver for the NACA0020 wing (same workflow as
FinRay/run_case.py). For each case, common/run_records.run_case():

  1. writes case_params.json and patches the case's "input2d" values into
     input2d
  2. runs Wing_Geom.py, then main2d.py, each in a fresh subprocess
  3. moves viz_IB2d/ (with the hooks' wing_history.npz / fluid_history.npz)
     into a new folder results/<case_name>_<YYYYmmdd-HHMMSS>/, with copies of
     case_params.json, input2d, wing.vertex, wing.target and wing_mesh.npz
  4. writes that folder's run_info.json (inputs, mesh metrics, timing,
     outcome, code version), then runs postprocess.py on it (the comparison
     with the SoRoLily reference: analysis section, comparison.png, ...)
  5. restores case_params.json and input2d

The sweeps are defined in sweeps.py and chosen with --sweep. On the cluster,
slurm/submit_sweep.sh runs one case per Slurm array task (--list, --case).

Run from this folder:
    python run_case.py [--sweep NAME] [--list | --case N]
    python run_case.py index [--all] [DIR ...]      (rebuild the index only)
'''

import json
import signal
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
import run_records as rr
import sweeps

EXAMPLE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = EXAMPLE_DIR / "results"

# Columns of index.csv (every field is in run_info.json and index_full.csv:
# `python run_case.py index --all`).
INDEX_COLUMNS = [
    # identity
    "folder", "run_id", "inputs.sweep", "timing.start", "backfilled",
    "provenance.commit", "provenance.branch", "provenance.uncommitted_changes",
    # parameters
    "inputs.geom.fill_density_multiplier", "inputs.geom.stiffness_scale",
    "inputs.geom.target_stiffness_factor", "inputs.geom.E_flex", "inputs.geom.init_q",
    "inputs.flow.U", "inputs.flow.fringe_rate", "inputs.flow.mean_rate", "inputs.flow.start_ramp",
    "inputs.input2d.mu", "inputs.input2d.rho", "inputs.input2d.dt", "inputs.input2d.Tfinal",
    "inputs.input2d.Nx", "inputs.input2d.Lx", "inputs.input2d.print_dump",
    # mesh
    "mesh.edge_spacing", "mesh.n_points", "mesh.n_springs", "mesh.stiffness_scale",
    "mesh.stiffness_scale_measured", "mesh.k_max", "mesh.n_target_points", "mesh.k_Target",
    # timing
    "timing.sim_wall_time_s", "timing.sim_time_per_step_ms",
    # outcome
    "outcome.status", "outcome.first_nan_time", "outcome.last_output_time", "outcome.interrupted",
    # analysis (postprocess.py)
    "analysis.stable", "analysis.first_nonfinite_time", "analysis.t_end", "analysis.re_chord", "analysis.clamp_drift_max", "analysis.force_ringing_rms",
    "analysis.te.y_initial", "analysis.te.t50", "analysis.te.t10", "analysis.te.y_min_release",
    "analysis.te.y_mean_0deg", "analysis.te.y_std_0deg", "analysis.te.y_mean_15deg",
    "analysis.te.y_std_15deg", "analysis.te.rms_diff_release", "analysis.te.rms_diff_0deg",
    "analysis.te.rms_diff_15deg",
    "analysis.force.fx_mean_0deg", "analysis.force.fy_mean_0deg",
    "analysis.force.fx_mean_15deg", "analysis.force.fy_mean_15deg",
    "analysis.force.fy_std_15deg", "analysis.force.rms_diff_fy_15deg",
    "analysis.shape.rms_dist_t0p25", "analysis.shape.rms_dist_t0p5", "analysis.shape.rms_dist_t1",
    "analysis.shape.rms_dist_mean_0deg", "analysis.shape.rms_dist_mean_15deg",
    "analysis.pressure.rms_diff_mean_15deg", "analysis.pressure.n_frames",
    "analysis.fluid.u_mean_end", "analysis.fluid.gs_ke_growth_rate",
    "analysis.fluid.gs_ke_fraction_end", "analysis.fluid.u_max", "analysis.fluid.cell_re_max",
]

BASE_CASE = json.loads((EXAMPLE_DIR / "case_params.json").read_text())
POST_SCRIPTS = ("postprocess.py",)
COPY_FILES = ("wing.vertex", "wing.target", "wing_mesh.npz")


def _interrupt(signum, frame):
    raise KeyboardInterrupt


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--sweep", default="base", metavar="NAME",
                        help=f"sweep to run, from sweeps.py: {', '.join(sweeps.SWEEPS)}")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--list", action="store_true", help="print the number and name of each case, then exit")
    group.add_argument("--case", type=int, metavar="N", help="run only case N of the sweep (0-based)")
    args = parser.parse_args()

    cases = sweeps.build(args.sweep, BASE_CASE)
    if args.list:
        for i, case in enumerate(cases):
            print(i, case["case_name"])
        return
    if args.case is not None:
        if not 0 <= args.case < len(cases):
            sys.exit(f"--case {args.case} out of range: the sweep has {len(cases)} case(s)")
        cases = [cases[args.case]]

    signal.signal(signal.SIGTERM, _interrupt)  # Slurm / an IDE's stop button: treat like Ctrl+C
    try:
        for case in cases:
            rr.run_case(case, EXAMPLE_DIR, RESULTS_DIR, "Wing_Geom.py",
                        post_scripts=POST_SCRIPTS, copy_files=COPY_FILES)
    finally:
        if args.case is None and RESULTS_DIR.exists():
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
