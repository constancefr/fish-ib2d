'''
Run records for IB2d cases. Every run folder gets a run_info.json (inputs,
mesh, timing, outcome, code provenance) -- the source of truth for that run --
and an index is generated from those files, one row per run, including runs
moved into subfolders (e.g. results/wrong_actuation/...; the `folder` column
gives each run's path): <results dir>/index.csv with a case's chosen columns
(see INDEX_COLUMNS in cases/FinRay/run_case.py), or index_full.csv with every
field. Indexes are never edited by hand: rebuild them whenever runs are added
or moved.

run_case() runs one case of a sweep (geometry script, then main2d.py, each
in a fresh subprocess), archives it and writes its record; each case's
run_case.py only defines its sweep and its INDEX_COLUMNS.

CLI (paths relative to where it is run from):
    python run_records.py index RESULTS_DIR [...]   # rebuild RESULTS_DIR/index_full.csv (every field)
    python run_records.py backfill RUN_DIR [...]    # records for runs made before records existed
'''

import csv
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime
from importlib import metadata
from pathlib import Path

PACKAGES = ('numpy', 'scipy', 'matplotlib', 'meshpy', 'numba')
MAX_UNTRACKED_BYTES = 200_000  # untracked files larger than this are listed, not copied into code.diff
# generated from the code + case_params.json (both recorded), so left out of code.diff
GENERATED = ('*.vertex', '*.spring', '*.target', '*.geo_connect', '*.beam', '*.png')


def read_input2d(path):
    '''Every `name = value` entry of an input2d file, numbers as int/float.'''
    out = {}
    for m in re.finditer(r'^\s*(\w+)\s*=\s*([^\s%]+)', Path(path).read_text(), re.MULTILINE):
        value = m.group(2).strip('"')
        for cast in (int, float):
            try:
                value = cast(value)
                break
            except ValueError:
                pass
        out[m.group(1)] = value
    return out


def _git(*args, cwd):
    return subprocess.run(['git', *args], cwd=cwd, capture_output=True, text=True).stdout


def provenance(case_dir, diff_path):
    '''
    The code a run used: commit, branch, submodule commits, package versions,
    machine. Uncommitted changes (tracked diffs, plus small untracked files)
    are saved to diff_path, so a run from a dirty tree is still reproducible.
    Generated files (GENERATED) are left out of the diff.
    '''
    status = _git('status', '--porcelain', cwd=case_dir).splitlines()
    if status:
        root = Path(_git('rev-parse', '--show-toplevel', cwd=case_dir).strip())
        diff = _git('diff', 'HEAD', '--', '.', *(f':!{g}' for g in GENERATED), cwd=root)
        for f in _git('ls-files', '--others', '--exclude-standard', cwd=root).splitlines():
            if any(Path(f).match(g) for g in GENERATED):
                continue
            if (root / f).stat().st_size <= MAX_UNTRACKED_BYTES:
                diff += _git('diff', '--no-index', '/dev/null', f, cwd=root)
        Path(diff_path).write_text(diff)

    versions = {}
    for p in PACKAGES:
        try:
            versions[p] = metadata.version(p)
        except metadata.PackageNotFoundError:
            versions[p] = None
    return dict(
        commit=_git('rev-parse', 'HEAD', cwd=case_dir).strip(),
        branch=_git('rev-parse', '--abbrev-ref', 'HEAD', cwd=case_dir).strip(),
        uncommitted_changes=bool(status),
        changed_files=[s.strip() for s in status],
        diff_file=Path(diff_path).name if status else None,
        submodules=_git('submodule', 'status', cwd=case_dir).strip(),
        python=platform.python_version(),
        packages=versions,
        machine=platform.node(),
        platform=platform.platform(),
    )


def outcome(viz_dir, dt, print_dump, t_final):
    '''
    How far a run got, from its lagPtsConnect.*.vtk frames (frame k is
    written at t = k * print_dump * dt). NaN never clears once it appears,
    so the first NaN frame is found by bisection.
    '''
    frames = sorted(Path(viz_dir).glob('lagPtsConnect.*.vtk'))
    if not frames:
        return dict(status='no output', first_nan_time=None, last_output_time=None, n_frames=0)

    def has_nan(f):
        return 'nan' in f.read_text().lower()

    def time_of(f):
        return round(int(f.name.split('.')[1]) * print_dump * dt, 9)

    first_nan = None
    if has_nan(frames[-1]):
        lo, hi = -1, len(frames) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            lo, hi = (lo, mid) if has_nan(frames[mid]) else (mid, hi)
        first_nan = time_of(frames[hi])

    last = time_of(frames[-1])
    if first_nan is not None:
        status = 'diverged'
    elif last >= t_final - print_dump * dt - 1e-9:
        status = 'completed'
    else:
        status = 'incomplete'
    return dict(status=status, first_nan_time=first_nan, last_output_time=last, n_frames=len(frames))


def inputs(case_params, input2d):
    '''The run's inputs: case_params sections (geom, actuation, ...) plus every
    input2d value IB2d read (which supersedes case_params["input2d"]).'''
    return {**{k: v for k, v in case_params.items() if k not in ('case_name', 'input2d')},
            'input2d': input2d}


def write_run_info(run_dir, info):
    (Path(run_dir) / 'run_info.json').write_text(json.dumps(info, indent=2) + '\n')


def _flatten(d, prefix=''):
    out = {}
    for k, v in d.items():
        if isinstance(v, dict):
            out.update(_flatten(v, f'{prefix}{k}.'))
        elif isinstance(v, list):
            out[f'{prefix}{k}'] = '; '.join(map(str, v))
        else:
            out[f'{prefix}{k}'] = v
    return out


def write_index(results_dir, columns=None):
    '''
    One row per run_info.json at any depth below results_dir, plus `folder`
    (the run's path relative to results_dir, which shows any grouping into
    subfolders). With `columns` (flattened keys such as "inputs.input2d.dt"),
    writes results_dir/index.csv with just those, in that order; without,
    results_dir/index_full.csv with every field.
    '''
    results_dir = Path(results_dir)
    rows = []
    for root, dirs, files in os.walk(results_dir):
        if 'run_info.json' in files:
            dirs.clear()  # a run folder: don't descend into its output
            rows.append({'folder': str(Path(root).relative_to(results_dir)),
                         **_flatten(json.loads((Path(root) / 'run_info.json').read_text()))})
    rows.sort(key=lambda r: str(r.get('timing.start') or ''))
    if columns is None:
        columns, path = list(dict.fromkeys(k for r in rows for k in r)), results_dir / 'index_full.csv'
    else:
        path = results_dir / 'index.csv'
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, columns, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    return path


def backfill(run_dir):
    '''
    run_info.json for a run made before records existed: inputs from its
    saved case_params.json / input2d, outcome from its output. Mesh and code
    version are unknown. Start time and sim wall time are estimated from the
    first and last output files' modification times, which is only valid if
    the folder was moved, not copied (an output span under 1 s means the
    times were reset by a copy, so they are dropped). Folders without an input2d are
    skipped (returns None).
    '''
    run_dir = Path(run_dir)
    if not (run_dir / 'input2d').exists():
        return None
    params_path = run_dir / 'case_params.json'
    params = json.loads(params_path.read_text()) if params_path.exists() else {}
    inp = read_input2d(run_dir / 'input2d')
    frames = sorted((run_dir / 'viz_IB2d').glob('lagPtsConnect.*.vtk'))
    t0, t1 = (frames[0].stat().st_mtime, frames[-1].stat().st_mtime) if frames else (None, None)
    if t0 is not None and t1 - t0 < 1.0:
        t0 = t1 = None
    info = dict(
        run_id=run_dir.name,
        case_name=params.get('case_name', run_dir.name),
        backfilled=True,
        inputs=inputs(params, inp),
        mesh=None,
        timing=dict(start=datetime.fromtimestamp(t0).isoformat(timespec='seconds') if t0 else None,
                    sim_wall_time_s_approx=round(t1 - t0, 1) if t0 else None),
        outcome=outcome(run_dir / 'viz_IB2d', inp['dt'], inp['print_dump'], inp['Tfinal']),
        provenance=None,
    )
    write_run_info(run_dir, info)
    return info



# --- Running a case -------------------------------------------------------

def patch_input2d(path, values):
    '''Rewrites matching `name = value` lines of the input2d file at `path`
    in place, keeping comments and ordering.'''
    path = Path(path)
    text = path.read_text()
    for key, val in values.items():
        # Matches "key = 0.01" or "key= 0.01   % comment"; keeps the comment.
        pattern = re.compile(rf"^(\s*{re.escape(key)}\s*=\s*)([^\s%]+)", re.MULTILINE)
        if not pattern.search(text):
            raise ValueError(
                f"Couldn't find '{key}' in input2d -- check the exact "
                f"parameter name/formatting in your input2d file.")
        text = pattern.sub(lambda m: m.group(1) + repr(val), text)
    path.write_text(text)


def run_case(case, example_dir, results_dir, geom_script, post_scripts=(), copy_files=()):
    '''
    Runs one case into its own new folder results_dir/<case_name>_<YYYYmmdd-HHMMSS>/
    (never overwrites an earlier run) and writes its run_info.json there:
    inputs, mesh metrics, timing, outcome, and the code version (commit,
    plus code.diff if the tree had uncommitted changes).

      1. writes case_params.json (read by the case's geometry and runtime
         scripts) and patches case["input2d"] into input2d
      2. runs geom_script, then main2d.py, each in a fresh subprocess (a
         case's runtime scripts may cache state in module globals, e.g.
         FinRay's update_Target_Point_Positions.py caches its baseline, so
         every case needs a clean interpreter, exactly as when running
         `python main2d.py` by hand). geom_script must write mesh_info.json.
      3. moves viz_IB2d/ and hier_IB2d_data/ into the run folder, with
         copies of case_params.json, input2d and copy_files
      4. writes run_info.json, then runs each post_script as
         `python <post_script> <run folder>` (a failure is reported, not
         raised)

    A run that is stopped (Ctrl+C, or an IDE's stop button: the caller should
    turn SIGTERM into KeyboardInterrupt) is still archived and recorded, with
    outcome.interrupted = True, and then KeyboardInterrupt is re-raised to
    stop the sweep. Returns the run folder.
    '''
    example_dir, results_dir = Path(example_dir), Path(results_dir)
    for folder in ("viz_IB2d", "hier_IB2d_data"):
        if (example_dir / folder).exists():
            raise RuntimeError(
                f"{folder}/ is left over from an earlier run that stopped before it was "
                f"archived. Move it into that run's results folder (then `python "
                f"../../common/run_records.py backfill <that folder>`), or delete it.")
    start = datetime.now()
    run_id = f"{case['case_name']}_{start:%Y%m%d-%H%M%S}"
    case_dir = results_dir / run_id
    case_dir.mkdir(parents=True)
    print(f"\n=== Running case: {run_id} ===")
    # code version as it is at launch (edits made while a sweep runs don't count)
    code = provenance(example_dir, case_dir / "code.diff")

    params_path = example_dir / "case_params.json"
    params_path.write_text(json.dumps(case, indent=2))
    patch_input2d(example_dir / "input2d", case["input2d"])

    env = dict(os.environ, MPLBACKEND="Agg")  # plt.show() is a no-op and doesn't block
    t0 = time.perf_counter()
    subprocess.run([sys.executable, geom_script], cwd=example_dir, env=env, check=True)
    t_geom = time.perf_counter() - t0

    # If interrupted, subprocess.run stops main2d.py; the run is still archived below.
    t0 = time.perf_counter()
    try:
        exit_code = subprocess.run([sys.executable, "main2d.py"], cwd=example_dir).returncode
        interrupted = False
    except KeyboardInterrupt:
        exit_code, interrupted = None, True
    t_sim = time.perf_counter() - t0

    for folder in ("viz_IB2d", "hier_IB2d_data"):
        src = example_dir / folder
        if src.exists():
            shutil.move(str(src), str(case_dir / folder))
    for name in ("case_params.json", "input2d", *copy_files):
        shutil.copy(example_dir / name, case_dir / name)

    inp = read_input2d(case_dir / "input2d")
    n_steps = round(inp["Tfinal"] / inp["dt"])
    info = dict(
        run_id=run_id,
        case_name=case["case_name"],
        inputs=inputs(case, inp),
        mesh=json.loads((example_dir / "mesh_info.json").read_text()),
        timing=dict(start=start.isoformat(timespec="seconds"),
                    end=datetime.now().isoformat(timespec="seconds"),
                    geometry_wall_time_s=round(t_geom, 1),
                    sim_wall_time_s=round(t_sim, 1),
                    n_steps=n_steps,
                    sim_time_per_step_ms=None if interrupted else round(1e3 * t_sim / n_steps, 2)),
        outcome=dict(outcome(case_dir / "viz_IB2d", inp["dt"], inp["print_dump"], inp["Tfinal"]),
                     exit_code=exit_code, interrupted=interrupted),
        provenance=code,
    )
    write_run_info(case_dir, info)

    status = info["outcome"]["status"]
    if status == "diverged":
        print(f"=== WARNING: {run_id} diverged (NaN from t = "
              f"{info['outcome']['first_nan_time']:.2f} s) -- see {case_dir} ===")
    else:
        print(f"=== {status}: results in {case_dir} ===")
    if interrupted:
        raise KeyboardInterrupt  # stop the sweep (the caller still rebuilds the index)
    if exit_code:
        raise subprocess.CalledProcessError(exit_code, "main2d.py")
    for script in post_scripts:
        r = subprocess.run([sys.executable, script, str(case_dir)], cwd=example_dir, env=env)
        if r.returncode:
            print(f"=== WARNING: {script} failed (exit {r.returncode}) for {run_id} ===")
    return case_dir

if __name__ == '__main__':
    if len(sys.argv) < 3 or sys.argv[1] not in ('index', 'backfill'):
        sys.exit(__doc__)
    for p in sys.argv[2:]:
        if sys.argv[1] == 'index':
            print('wrote', write_index(p))
        else:
            info = backfill(p)
            print(p, '->', info['outcome']['status'] if info else 'skipped (no input2d)')
