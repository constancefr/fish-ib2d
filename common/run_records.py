'''
Run records for IB2d cases. Every run folder gets a run_info.json (inputs,
mesh, timing, outcome, code provenance) -- the source of truth for that run --
and <results dir>/index.csv is generated from those files, one row per run.
The index is never edited by hand: rebuild it whenever runs are added.

CLI (paths relative to where it is run from):
    python run_records.py index RESULTS_DIR [...]   # rebuild RESULTS_DIR/index.csv
    python run_records.py backfill RUN_DIR [...]    # records for runs made before records existed
'''

import csv
import json
import platform
import re
import subprocess
import sys
from datetime import datetime
from importlib import metadata
from pathlib import Path

PACKAGES = ('numpy', 'scipy', 'matplotlib', 'meshpy', 'numba')
MAX_UNTRACKED_BYTES = 200_000  # untracked files larger than this are listed, not copied into code.diff


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
    '''
    status = _git('status', '--porcelain', cwd=case_dir).splitlines()
    if status:
        diff = _git('diff', 'HEAD', cwd=case_dir)
        root = Path(_git('rev-parse', '--show-toplevel', cwd=case_dir).strip())
        for f in _git('ls-files', '--others', '--exclude-standard', cwd=root).splitlines():
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


def write_index(results_dir):
    '''results_dir/index.csv: one row per run_info.json below it, every field a column.'''
    results_dir = Path(results_dir)
    rows = [_flatten(json.loads(p.read_text())) for p in sorted(results_dir.glob('*/run_info.json'))]
    rows.sort(key=lambda r: str(r.get('timing.start') or ''))
    columns = list(dict.fromkeys(k for r in rows for k in r))
    path = results_dir / 'index.csv'
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, columns)
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


if __name__ == '__main__':
    if len(sys.argv) < 3 or sys.argv[1] not in ('index', 'backfill'):
        sys.exit(__doc__)
    for p in sys.argv[2:]:
        if sys.argv[1] == 'index':
            print('wrote', write_index(p))
        else:
            info = backfill(p)
            print(p, '->', info['outcome']['status'] if info else 'skipped (no input2d)')
