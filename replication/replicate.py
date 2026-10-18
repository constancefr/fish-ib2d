'''
Re-runs the simulations listed in manifest.json, each with the exact code
that first produced it (see README.md).

For each run: exports the run's commit into a private source tree (git
archive, plus any files that were untracked at the time), links in this
checkout's external/IB2d (same pinned commit for every run), copies the run's
own case_params.json and input2d into the case folder, then runs the geometry
script and main2d.py in fresh subprocesses, as run_case.py does. Output goes
to <out>/<run path>/ (the original layout, e.g.
results_replicated/FinRay/results_approach2/mult32.00_scale1.00_20261014-114157/)
with a run_info.json in the usual format plus a "replication" section that
compares the new run with the original (initial points, positions at a
reference frame, first NaN time).

Usage (from the repo root):
    python replication/replicate.py list [--group G ...]
    python replication/replicate.py run RUN_PATH ...            # paths as listed
    python replication/replicate.py run --group approach1       # a whole group
    python replication/replicate.py run --all                   # everything, in order
    python replication/replicate.py run --all --index 7         # 8th run of the selection (job arrays)
Options for run:
    --out DIR      output root (default: results_replicated/ in the repo root)
    --geom-only    build the geometry and check it against the original; no simulation
    --force        rerun a run whose output folder already exists
    --keep-src     keep the exported source tree (<out>/_work/<run path>/src)
'''

import argparse
import io
import json
import os
import platform
import shutil
import signal
import subprocess
import sys
import tarfile
import time
from datetime import datetime
from importlib import metadata
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(REPO / "common"))
import run_records as rr  # noqa: E402

MANIFEST = json.loads((HERE / "manifest.json").read_text())
RUNS = MANIFEST["runs"]
PACKAGES = (*rr.PACKAGES, "vtk", "pyFFTW")  # IB2d uses pyFFTW for the fluid solve when installed


def _git(*args, **kw):
    return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, check=True, **kw).stdout


def vtk_points(path):
    '''x, y of every point in a legacy ASCII .vtk file (IB2d writes 2.0 or 5.1).'''
    lines = Path(path).read_text().split("\n")
    i = next(k for k, line in enumerate(lines) if line.startswith("POINTS"))
    n = int(lines[i].split()[1])
    vals, k = [], i + 1
    while len(vals) < 3 * n:
        vals += lines[k].split()
        k += 1
    return np.array(vals[:3 * n], float).reshape(n, 3)[:, :2]


def vertex_points(path):
    return np.loadtxt(path, skiprows=1).reshape(-1, 2)


def max_diff(a, b):
    if a.shape != b.shape:
        return f"shape {a.shape} vs {b.shape}"
    return float(np.abs(a - b).max())


def export_source(entry, src):
    '''The run's commit as a plain tree at src, with external/IB2d linked in.'''
    tar = _git("archive", "--format=tar", entry["commit"])
    with tarfile.open(fileobj=io.BytesIO(tar)) as t:
        t.extractall(src, filter="data")
    for extra in entry["extra_files"]:
        (src / extra["path"]).parent.mkdir(parents=True, exist_ok=True)
        (src / extra["path"]).write_bytes(_git("show", f"{extra['from_commit']}:{extra['path']}"))

    pinned = _git("ls-tree", entry["commit"], "external/IB2d", text=True).split()[2]
    have = _git("-C", "external/IB2d", "rev-parse", "HEAD", text=True).strip()
    if have != pinned:
        sys.exit(f"external/IB2d is at {have}, but {entry['commit'][:7]} pins {pinned}: "
                 f"run `git submodule update --init` first")
    ib2d = src / "external" / "IB2d"
    if ib2d.exists():
        shutil.rmtree(ib2d)
    ib2d.parent.mkdir(parents=True, exist_ok=True)
    ib2d.symlink_to(REPO / "external" / "IB2d", target_is_directory=True)


def provenance(entry):
    versions = {}
    for p in PACKAGES:
        try:
            versions[p] = metadata.version(p)
        except metadata.PackageNotFoundError:
            versions[p] = None
    return dict(
        commit=entry["commit"],
        branch=None,
        uncommitted_changes=bool(entry["extra_files"]),
        changed_files=[e["path"] for e in entry["extra_files"]],
        diff_file=None,
        submodules=f"{MANIFEST['ib2d_submodule']} external/IB2d",
        replication_kit_commit=_git("rev-parse", "HEAD", text=True).strip(),
        python=platform.python_version(),
        packages=versions,
        machine=platform.node(),
        platform=platform.platform(),
        slurm_job=os.environ.get("SLURM_JOB_ID"),
    )


def run_one(entry, out_root, geom_only=False, force=False, keep_src=False):
    run_dir = out_root / entry["path"]
    done = run_dir / "run_info.json"
    if done.exists() and not force and not geom_only:
        if not json.loads(done.read_text())["outcome"]["interrupted"]:
            print(f"--- skip (already done): {entry['path']}")
            return
        print(f"--- redoing (stopped early last time, e.g. a job time limit): {entry['path']}")
        force = True
    if run_dir.exists() and not geom_only:
        if not force:
            sys.exit(f"{run_dir} exists without a run_info.json (an earlier attempt stopped "
                     f"early?): delete it, or pass --force")
        shutil.rmtree(run_dir)

    work = out_root / "_work" / entry["path"]
    if work.exists():
        shutil.rmtree(work)
    src = work / "src"
    src.mkdir(parents=True)
    print(f"\n=== {entry['path']}  ({entry['group']}, code {entry['commit'][:7]}) ===")
    export_source(entry, src)

    inputs_dir = HERE / "runs" / entry["path"]
    case_dir = src / "cases" / entry["case"]
    for name in ("case_params.json", "input2d"):
        shutil.copy(inputs_dir / name, case_dir / name)
    ref = np.load(inputs_dir / "reference_points.npz")
    inp = rr.read_input2d(case_dir / "input2d")
    struct = inp["string_name"]

    env = dict(os.environ, MPLBACKEND="Agg")
    start = datetime.now()
    t0 = time.perf_counter()
    subprocess.run([sys.executable, entry["geom_script"]], cwd=case_dir, env=env, check=True,
                   stdout=subprocess.DEVNULL if geom_only else None)
    t_geom = time.perf_counter() - t0
    d0 = max_diff(vertex_points(case_dir / f"{struct}.vertex"), ref["frame0"])
    print(f"initial points vs original: max diff {d0}")
    if geom_only:
        if not keep_src:
            shutil.rmtree(work)
        return d0

    t0 = time.perf_counter()
    try:
        exit_code, interrupted = subprocess.run([sys.executable, "main2d.py"], cwd=case_dir).returncode, False
    except KeyboardInterrupt:
        exit_code, interrupted = None, True
    t_sim = time.perf_counter() - t0

    run_dir.mkdir(parents=True)
    for folder in ("viz_IB2d", "hier_IB2d_data"):
        if (case_dir / folder).exists():
            shutil.move(str(case_dir / folder), str(run_dir / folder))
    for name in ("case_params.json", "input2d", "mesh_info.json", *entry["copy_files"]):
        if (case_dir / name).exists():
            shutil.copy(case_dir / name, run_dir / name)

    out = rr.outcome(run_dir / "viz_IB2d", inp["dt"], inp["print_dump"], inp["Tfinal"])
    ref_file = run_dir / "viz_IB2d" / f"lagsPts.{int(ref['ref_index']):04d}.vtk"
    orig = entry["original"]
    check = dict(
        initial_points_max_diff=d0,
        ref_frame=int(ref["ref_index"]), ref_time=float(ref["ref_time"]),
        ref_points_max_diff=max_diff(vtk_points(ref_file), ref["ref"]) if ref_file.exists() else "frame missing",
        status=[orig["status"], out["status"]],
        first_nan_time=[orig["first_nan_time"], out["first_nan_time"]],
    )
    case_params = json.loads((run_dir / "case_params.json").read_text())
    n_steps = round(inp["Tfinal"] / inp["dt"])
    info = dict(
        run_id=Path(entry["path"]).name,
        case_name=case_params.get("case_name", Path(entry["path"]).name),
        inputs=rr.inputs(case_params, inp),
        mesh=json.loads((run_dir / "mesh_info.json").read_text()) if (run_dir / "mesh_info.json").exists() else None,
        timing=dict(start=start.isoformat(timespec="seconds"),
                    end=datetime.now().isoformat(timespec="seconds"),
                    geometry_wall_time_s=round(t_geom, 1),
                    sim_wall_time_s=round(t_sim, 1),
                    n_steps=n_steps,
                    sim_time_per_step_ms=None if interrupted else round(1e3 * t_sim / n_steps, 2)),
        outcome=dict(out, exit_code=exit_code, interrupted=interrupted),
        provenance=provenance(entry),
        replication=dict(of=entry["path"], group=entry["group"], original=orig, check=check),
    )
    rr.write_run_info(run_dir, info)
    print(f"=== {out['status']} (original: {orig['status']}); initial points diff {d0}, "
          f"frame {check['ref_frame']} diff {check['ref_points_max_diff']}; "
          f"first NaN {check['first_nan_time']} ===")

    if not interrupted and not exit_code:
        for script in entry["post_scripts"]:
            r = subprocess.run([sys.executable, script, str(run_dir)], cwd=case_dir, env=env)
            if r.returncode:
                print(f"=== WARNING: {script} failed (exit {r.returncode}) ===")
    if not keep_src:
        shutil.rmtree(work)
    if interrupted:
        raise KeyboardInterrupt
    if exit_code:
        raise subprocess.CalledProcessError(exit_code, "main2d.py")


def select(args):
    if args.all:
        chosen = list(RUNS)
    elif args.group:
        chosen = [e for e in RUNS if e["group"] in args.group]
    else:
        by_path = {e["path"]: e for e in RUNS}
        unknown = [p for p in args.runs if p.rstrip("/") not in by_path]
        if unknown:
            sys.exit(f"not in manifest.json: {unknown} (see `replicate.py list`)")
        chosen = [by_path[p.rstrip("/")] for p in args.runs]
    if not chosen:
        sys.exit("no runs selected")
    if args.index is not None:
        chosen = [chosen[args.index]]
    return chosen


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    ls = sub.add_parser("list")
    ls.add_argument("--group", nargs="+")
    r = sub.add_parser("run")
    r.add_argument("runs", nargs="*")
    r.add_argument("--all", action="store_true")
    r.add_argument("--group", nargs="+")
    r.add_argument("--index", type=int)
    r.add_argument("--out", type=Path, default=REPO / "results_replicated")
    r.add_argument("--geom-only", action="store_true")
    r.add_argument("--force", action="store_true")
    r.add_argument("--keep-src", action="store_true")
    args = p.parse_args()

    if args.cmd == "list":
        runs = [e for e in RUNS if not args.group or e["group"] in args.group]
        for i, e in enumerate(runs):
            q, o = e["params"], e["original"]
            print(f"{i:3d}  {e['group']:27s} {e['commit'][:7]}  Nx={q['Nx']:<3} dt={q['dt']:<7g} "
                  f"T={q['Tfinal']:<4g} {o['status']:9s} {o['sim_wall_time_s'] or '?':>7}s  {e['path']}")
        print(f"{len(runs)} runs")
        return

    signal.signal(signal.SIGTERM, lambda *_: (_ for _ in ()).throw(KeyboardInterrupt()))
    args.out = args.out.resolve()
    failed = []
    for entry in select(args):
        try:
            run_one(entry, args.out, args.geom_only, args.force, args.keep_src)
        except subprocess.CalledProcessError as e:
            print(f"!!! FAILED: {entry['path']}: {e}")
            failed.append(entry["path"])
    work = args.out / "_work"
    if work.exists() and not any(p.is_file() for p in work.rglob("*")):
        shutil.rmtree(work)  # only empty folders left
    if failed:
        sys.exit(f"{len(failed)} failed: {failed}")


if __name__ == "__main__":
    main()
