# FinRay: experiment plan (numerical verification, then physics)

Drafted 2026-10-18 from a planning discussion; to be carried out on the cluster
(`slurm/submit_sweep.sh`, 8 parallel jobs). Read `NOTES.md` first: it holds the
model, the pipeline and every finding this plan builds on.

**Ground rule: the user decides every aspect of the experimental setup.**
Items marked **OPEN** are not decided. The defaults given are suggestions only:
ask the user before choosing, and record each decision in the "Decisions" table.

## 1. Goal and approach

Find numerical settings where the results no longer depend on the numerics
("verification"), at the lowest cost, then study the physical parameters at
those settings. Numerical parameters are refined one at a time (convergence is
one-dimensional per parameter), except the three pairs expected to interact:
density × dt (stiffer meshes need smaller dt), `Nx` × `mu` (both set the cell
Reynolds number), and density × `Nx` (mesh points per grid cell).

| kind | parameters | aim |
|---|---|---|
| numerical (should not change the answer) | density (`ribbon_fill_density_multiplier`), dt, `Nx`, `target_stiffness_factor`, `STIFFNESS_RAMP_TIME` | show convergence, then pick the cheapest sufficient settings |
| physical (the study) | `FREQUENCY`, `AMPLITUDE`, `E_material`, `mu` (if treated as a model choice) | explore their effects at verified settings |

In IB2d, ds = Lx/(2 Nx) is fixed by the grid, so "grid resolution" and "ds" are
one parameter. The free one is the mesh spacing relative to the grid cell (h/dx).

## 2. Decisions so far (user, 2026-10-18)

| # | decision |
|---|---|
| 1 | **Regime: the real actuator's.** Currently `FREQUENCY` = 0.5 Hz (one side pushes and pulls back in 1 s, then the other side in the next 1 s: a 2 s cycle, which is exactly what `update_Target_Point_Positions.py` does at 0.5 Hz) and `AMPLITUDE` 2.5-3 cm. These may change with hardware experiments, so **the whole plan must be re-runnable from config** (one place to change frequency and amplitude, run lengths derived from the frequency). |
| 2 | **Convergence is judged on the tip trajectory.** |
| 3 | **Viscosity: as close to water as possible.** A viscosity sweep is fine, but the aim is to stay near water, not to pick an inflated `mu` for convenience. |
| 4 | **Compute: the cluster, 8 jobs in parallel** (one CPU each; IB2d is single-threaded). |
| 5 | Stages 2 and 4 of the first draft (stability, then coupled grid convergence) are **merged into one grid ladder** (stage C below): each run gives both the stability measure and the tip trajectory. |

## 3. Facts the plan relies on (details in NOTES.md)

- **Springs are ~3× too soft at `ribbon_fill_stiffness_scale` = 1**: tributary-area springs give an effective modulus ~E/3 (NOTES §4 issue 1; the wing case calibrated scale ≈ 4 against a reference deflection). Every result depends on this, so calibration comes first (stage B).
- **dt limit at Nx = 64** (1 s runs, old regime 0.1 Hz / 1 cm): ×64 and ×128 stable at dt 4e-4, NaN at 5e-4 (×254 NaN at 5e-4, not bisected). At Nx = 32, ×256-equivalent meshes run at 5e-4. The failures come as the target stiffness finishes ramping in (0.2 s). NOTES §4 issue 2.
- **Fluid instability**: at water viscosity the grid-scale fluid kinetic energy grows steadily (central-difference advection at high cell Re). With the outer-face actuation patch, a full 15 s sweep at Nx = 32 finished at every density, but at ×128 the energy was running away by 15 s. The NACA case went NaN at cell Re ~1250 and was stable at 30× water viscosity. NOTES §4 issue 2.
- **Cell Re estimate at the new regime**: the actuator's peak speed is A·2π·f ≈ 0.03·π ≈ 0.09 m/s, so cell Re = u·dx/ν ≈ 2900 at Nx = 32, 1500 at Nx = 64, 730 at Nx = 128 (water). Expect water viscosity to need Nx ≥ 128, or `mu` raised a few times. The tip moves faster than the actuator, so these are lower bounds.
- **Cost per step** (M-series Mac, single core): Nx = 32 ×128 ≈ 14 ms; Nx = 64 ≈ 47 ms; Nx = 128 ≈ 4× Nx = 64 (estimate). IB2d's spring force is a Python loop (~6 µs per spring), so cost also grows with density. Measure on the cluster: every run records `timing.sim_time_per_step_ms`.
- **Mesh noise**: changing the density moves every Halton point, so results jitter on top of the convergence trend. Measure it with several Halton seeds at fixed density (needs a seed option, stage A).

## 4. Stages, in order

| stage | what | design | runs (est.) | output |
|---|---|---|---|---|
| **A. Tooling** | see list below | — | — | the measures and sweep definitions every later stage needs |
| **B. Structure only (no fluid): convergence + calibration** | density ×1 → ×1024 in factors of 4 (h halves each step); 3-5 Halton seeds each; static load | port the wing's fluid-free static solver (`cases/NACA0020_Wing/verify_structure.py`) | ~25 solves, seconds each | converged density; seed jitter; `ribbon_fill_stiffness_scale` against the strip test (§7: side wall + shortest ray, plane strain) |
| **C. Grid ladder: stability + grid convergence (merged)** | `Nx` ∈ {32, 64, 128} × `mu` ∈ {1, 3, 10, (30)} × water, at the stage-B density with h/dx fixed | runs long enough for the tip measures (**OPEN: number of cycles**, default 3 after the ramp ≈ 7-8 s at 0.5 Hz); dt from stage C0 | ~9-12, one or two batches of 8; Nx = 128 runs ~1-1.5 h each (estimate) | per (`Nx`, `mu`): stable or not (fluid energy growth rate), tip trajectory. Then: the `mu` closest to water that is stable at an affordable `Nx`, and the convergence order + error bar (GCI, 3 grids) of the tip trajectory at that `mu` |
| C0 (before C) | **critical dt** per `Nx` at the stage-B density (and one denser) at the new regime | bisection on a dt ladder whose 1/dt are multiples of `print_dump` (e.g. 4e-4, 3.125e-4, 2.5e-4, 2e-4, 1.25e-4, 1e-4); pass = no NaN within 1 s | ~3 per combination | dt_crit; production runs use a fixed fraction (**OPEN**, default 0.7) |
| **D. Mesh-noise error bar** | 3-5 Halton seeds at the chosen setting | — | 3-5 | spread of the tip measures; a refinement level only counts as converged if its change exceeds this spread |
| **E. Actuation numerics** | `target_stiffness_factor` (e.g. 3, 7, 15, 30) × `STIFFNESS_RAMP_TIME` | short runs | ~6-8 | patch tracking error vs dt_crit; pick the factor |
| **F. Physics** | `FREQUENCY` and `AMPLITUDE` one at a time around the baseline (~5 levels each), then a ~4×4 grid over the interesting region; possibly `E_material` (0.45 vs 0.74 MPa, NOTES §4 issue 1) | at the settings from B-E | ~10 + ~16 | the study (**OPEN: ranges, levels, whether E is included**) |

### Stage A: tooling to build first

**Status (2026-10-18): all five built, not yet run on the cluster.** Files:
`sweeps.py` (2, 3), `postprocess.py` (1), `halton_seed` in
`common/mesh_utils.py` / `FinRay_Geom.py` (4), `verify_structure.py` (5).
NOTES.md §2 describes them.

1. **Post-processing → `run_info.json`**: tip trajectory per cycle after the ramp (amplitude, phase lag relative to the actuation, peak-to-peak), the growth rate of grid-scale fluid kinetic energy (3×3 box-filter residual; the scan used for NOTES §4 issue 2), time of first NaN. Add these as index columns. A post-processing step that can be re-run on existing runs without re-simulating.
2. **General sweep definitions**: `make_sweep()` currently varies only density × stiffness scale (hard-coded lists). Each stage needs its own sweep (any `case_params` key, lists of values or a small grid), selectable when submitting, so a stage can be re-run unchanged after editing the regime.
3. **Regime in one place**: frequency and amplitude in `case_params.json["actuation"]` (already), `Tfinal` derived from the number of cycles (**OPEN**) and the frequency, so changing the hardware values reruns everything consistently.
4. **Halton seed option** in `mesh_utils.fill_polygon_with_halton` (scrambled Halton with a seed; the current fill is unscrambled), recorded in `run_info.json`.
5. **Fluid-free static solver for FinRay** (stage B), modelled on the wing's `verify_structure.py` (linear sparse static solve + IB2d force-routine cross-check).

## 5. Choosing the cheapest sufficient settings

Set the tolerance before running (**OPEN**, e.g. tip amplitude and phase within
2%). At each refinement step, compare with the finest level or a Richardson
extrapolation. The chosen setting is the coarsest one that:

1. has an error below the tolerance, and the error is larger than the stage-D seed spread (otherwise the seed spread is the real limit);
2. stays stable for the run length needed (fluid energy does not run away);
3. runs at dt ≤ the chosen fraction of dt_crit.

Cost is read from the index (`timing.sim_time_per_step_ms` × steps). Report the
convergence order and the grid convergence index (GCI, Roache) on the tip measures.

## 6. Open decisions (ask the user)

Still open (decided items are in §7):

| # | decision | default suggestion |
|---|---|---|
| O1 | which tip measure is judged (amplitude, phase lag, whole-trajectory RMS), and the tolerance. Needed before stage C | 2% on amplitude and phase |
| O3 | `E_material` for the study: 0.74 MPa (no known source) or 0.45 MPa (UCL database, Dragon Skin 30). Needed before C0 (it changes the stable dt) | 0.45 MPa |
| O4 | `mu` levels in stage C | 1, 3, 10 × water (30 × only if needed) |
| O6 | number of Halton seeds (stages B and D) | 3 (5 if the spread is close to the tolerance) |
| O7 | fraction of dt_crit for production runs | 0.7 |
| O8 | stage F ranges and levels (the baseline is decided) | f 0.25-1 Hz, A 1.5-3.5 cm, 5 levels |

## 7. Decisions log

Each decision with a plain explanation, for readers new to simulation.

| date | # | decision | what it means, and why |
|---|---|---|---|
| 2026-10-18 | 1-5 | see §2 | planning discussion |
| 2026-10-18 | O9 | **Sweeps are Python functions** in `sweeps.py`, one per sweep, chosen by name (`run_case.py --sweep NAME`, `submit_sweep.sh <folder> NAME`). | A "sweep" is a batch of runs that differ in a few settings. Code (rather than a data file) makes linked settings easy, e.g. a smaller time step for each finer grid. |
| 2026-10-18 | O5 | **Run length = 1 start-up cycle + 3 measured cycles**: `Tfinal` = 4 / frequency (8 s at 0.5 Hz), worked out automatically unless a sweep sets `Tfinal` itself. | One cycle = both actuators push and return. The first cycle is thrown away because the fin and the water start from rest (and the actuator stiffness ramps in), so it isn't yet repeating. 3 cycles show whether the motion repeats. |
| 2026-10-18 | O1 (part) | **The "tip" is the outer tip vertex**: the sharp outer corner of the fin's outline. All candidate measures are computed (§4, stage A item 1); which one is judged, and the tolerance, are still open. | The same physical point exists in every mesh, whatever the density or random seed, so runs can be compared directly. |
| 2026-10-18 | O8 (part) | **Baseline regime: 0.5 Hz, 3 cm** (`case_params.json`). | The real actuator's frequency, and the upper end of its 2.5-3 cm stroke. |
| 2026-10-18 | O2 | **Stiffness reference: the strip test.** A straight piece of wall (3.96 mm thick) is meshed exactly like the fin, clamped at one end and pushed at the other, without fluid (`verify_structure.py`). Its bending is compared with Timoshenko beam theory (the textbook formula, including shear). The ratio gives `ribbon_fill_stiffness_scale`. A whole-fin finite-element reference was not chosen (for now). | The spring mesh is about 3× too soft at scale 1 (NOTES §4 issue 1), so it must be calibrated against something known. A strip has an exact textbook answer; the wing case used the same method. |
| 2026-10-18 | B setup | **Two strips: the side wall (138.4 mm) and the shortest ray (16.9 mm)**, both centreline lengths. | Long strips mostly bend; short, stubby ones also shear (about 6% of the deflection here). If both give the same scale, one number fits the whole fin; if not, we learn by how much they differ. |
| 2026-10-18 | B setup | **Plane strain, ν = 0.499**: the reference stiffness is E/(1-ν²) ≈ 1.33 E. | The 2D model treats the fin as infinitely deep. A 27 mm deep, 4 mm thick wall bends like a plate: it can't contract sideways, which makes it stiffer than a thin beam ("plane stress"). ν (Poisson's ratio) = 0.499 is silicone's. |
| 2026-10-18 | base | **Base numerics in `case_params.json`: `Nx` 32, ×8, dt 5e-4** (placeholders). | `Nx` = grid cells per side, ×8 = mesh density, dt = time step. Each stage's sweep sets these explicitly; the base only serves the `base`/`smoke` test sweeps, so it is cheap and known stable (old regime). The previous base (`Nx` 64, ×254) went NaN at this dt. |
| 2026-10-18 | workflow | **The user submits cluster jobs.** Claude prepares the code and sweeps and hands over an ordered submission plan. | Keeps control of what reaches the shared cluster with the user. |
