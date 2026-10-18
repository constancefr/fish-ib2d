# Replication kit

Everything needed to re-run, on another machine, every simulation that existed locally on 2026-10-16 (32 runs, ~15 GB of output), each with **the exact code that first produced it**. Made before the local output was deleted.

| File | What it is |
|---|---|
| `manifest.json` | one entry per run: the commit to run it with, group, key parameters, the original outcome, and how the commit was identified |
| `runs/<path>/` | per run: the original `case_params.json` and `input2d` (used as-is), `original_run_info.json`, and `reference_points.npz` (the original Lagrangian points at frame 0 and at the last frame before any NaN) |
| `replicate.py` | runs any subset of the manifest |
| `replicate_array.sbatch` | SLURM job array, one task per run |
| `environment.txt` | package versions the originals used |
| `original_indexes/` | the original `index.csv` / `index_full.csv` files |
| `runs/legacy/` | inputs of the two old runs that can't be reproduced (below) |

## How it works

For each run, `replicate.py`:
1. exports the run's commit with `git archive` into `<out>/_work/<path>/src/` (one private tree per run, so array tasks never collide), adds any file that was untracked at the time (`extra_files`), and links in this checkout's `external/IB2d` (every run used the same pinned IB2d commit, a5136ac).
2. copies the run's own `case_params.json` and `input2d` into the case folder, so no regex patching is involved.
3. runs the geometry script, then `main2d.py`, each in a fresh subprocess, just like `run_case.py`.
4. moves `viz_IB2d/` and `hier_IB2d_data/` to `<out>/<path>/` (same layout as the originals, e.g. `results_replicated/FinRay/results_approach2/mult32.00_scale1.00_20261014-114157/`) and writes a normal `run_info.json`. It also adds a `replication` section that compares the new run with the original:
   - `initial_points_max_diff`: the generated structure against the original frame 0. Expect about 5e-12, the precision of the saved VTK files.
   - `ref_points_max_diff`: the positions at `ref_frame` (the last frame before any NaN) against the original.
   - `status` and `first_nan_time` as `[original, new]` pairs.

Indexes work as before, e.g. `cd cases/FinRay && python run_case.py index ../../results_replicated/FinRay/results_approach2`.

## On the cluster

```sh
git clone --recurse-submodules <repo URL> && cd fish-ib2d
conda create -n ib2d python=3.14 && conda activate ib2d
pip install -r replication/environment.txt          # meshpy needs a C++ compiler if no wheel exists
python replication/replicate.py list                 # 32 runs, numbered
python replication/replicate.py run --all --geom-only --out /tmp/geomcheck   # quick check (a few minutes): each line should say max diff ~5e-12
mkdir -p logs && sbatch replication/replicate_array.sbatch
```

The commits must be reachable from the clone: `main` contains every commit used here (`approach2` and `wing-approach2` were merged into it, then deleted). Set `OUT=/path/on/scratch` when submitting to write elsewhere. Then edit the environment block and the `#SBATCH` lines for your cluster. On an M-series Mac, original wall times ranged from 30 s to 32 min; `replicate.py list` shows each one.

Rows of `replicate.py list`, to submit one group with `sbatch --array=<rows>`:

| Rows | Group | Commit | Runs |
|---|---|---|---|
| 0-5 | `approach1` | 96af2b4 (×0, ×8, ×64, ×128) / 60ce2bf (×1, ×32) | per-ribbon mesh, old actuation patch, Nx 32, 15 s |
| 6-11 | `approach2_wrong_actuation` | 0e6666d | single mesh, old (whole-thickness) actuation patch, Nx 32, 15 s |
| 12-17 | `approach2_Nx32` | 906dce4 | single mesh, outer-face actuation, ×0-×128, Nx 32, 15 s |
| 18-19 | `approach2_Nx32_extra` | fb3147b | ×0 repeat and ×256, Nx 32, 15 s |
| 20-26 | `approach2_Nx64` | fb3147b | ×0-×128 and ×254, Nx 64, dt 5e-4, 15 s |
| 27-30 | `approach2_dt_critical_Nx64` | fb3147b | ×64 / ×128 at dt 2.5e-4 and 4e-4, Nx 64, 1 s |
| 31 | `wing` | 10e8e4d | NACA0020, Nx 192, 10 s (then `analyze_run.py`) |

## How each run's code was identified

- **Runs with a record** (`approach2*`, `wing`): commit + the run's saved `code.diff`. Applying each diff and comparing the geometry, actuation, solver and `common/` code showed the uncommitted changes never touched simulation code at run time. They only touched `run_case.py`, `case_params.json`, notes, record-keeping and generated files. The only exception is the afternoon of 2026-10-14, whose `FinRay_Geom.py` is byte-identical to the next commit, fb3147b. So each run maps to one clean commit.
- **approach1** (backfilled, no commit recorded): the geometry was regenerated at candidate commits and compared with each run's frame 0. ×0, ×8, ×64 and ×128 (run 2026-10-09) match 96af2b4 to 5e-12. That needs `common/ib_units.py`, which was untracked then and identical in every later commit, so it is taken from 60ce2bf. ×1 and ×32 (run 2026-10-12) match 60ce2bf, the cleaned-up approach1 that `approach1_trimmed` (f60d680) only renames. They don't match 96af2b4, because the cleanup changed the Halton points.
- **Check**: `replicate.py run --all --geom-only` reproduced every run's initial points to 5e-12. Two full runs (approach1 ×0, approach2 Nx32 ×0) were rerun on the original Mac. Both were **bit-identical** to the originals (position difference 0.0 at frames 1203 and 1499), and approach1 ×0 diverged at 12.04 s, as the original did.

## Expect small differences on the cluster

Same code and inputs, but a different CPU, BLAS, FFTW plan and compiler give floating-point differences from the first step. Stable runs should agree closely. Runs that blow up (most of `approach1`, `wrong_actuation`, and the Nx 64 sweep) amplify tiny differences, so their `first_nan_time` may shift somewhat. That is chaos, not a code mismatch: the `initial_points_max_diff` check confirms the code and inputs.

## Not reproducible

- `old_results/FinRay_results_0` (×0, ×1; ×8 only had movie frames): made around 2026-09-30 in the old IB2d fork from uncommitted code (its geometry imports `_mesh_utils`, which isn't in b656ebf). It is also the pre-units-fix model (rho = 1, springs 1280× too soft). Its inputs are in `runs/legacy/`.
- `old_results/NACA0200_results_0`: no inputs or record were saved. It came from the first, superseded NACA implementation (2026-10-08).
