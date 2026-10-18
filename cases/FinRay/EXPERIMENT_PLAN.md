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
| **B. Structure only (no fluid): convergence + calibration** | density ×1 → ×1024 in factors of 4 (h halves each step); 3-5 Halton seeds each; static load | port the wing's fluid-free static solver (`cases/NACA0020_Wing/verify_structure.py`) | ~25 solves, seconds each | converged density; seed jitter; `ribbon_fill_stiffness_scale` against a reference (**OPEN: which reference**) |
| **C. Grid ladder: stability + grid convergence (merged)** | `Nx` ∈ {32, 64, 128} × `mu` ∈ {1, 3, 10, (30)} × water, at the stage-B density with h/dx fixed | runs long enough for the tip measures (**OPEN: number of cycles**, default 3 after the ramp ≈ 7-8 s at 0.5 Hz); dt from stage C0 | ~9-12, one or two batches of 8; Nx = 128 runs ~1-1.5 h each (estimate) | per (`Nx`, `mu`): stable or not (fluid energy growth rate), tip trajectory. Then: the `mu` closest to water that is stable at an affordable `Nx`, and the convergence order + error bar (GCI, 3 grids) of the tip trajectory at that `mu` |
| C0 (before C) | **critical dt** per `Nx` at the stage-B density (and one denser) at the new regime | bisection on a dt ladder whose 1/dt are multiples of `print_dump` (e.g. 4e-4, 3.125e-4, 2.5e-4, 2e-4, 1.25e-4, 1e-4); pass = no NaN within 1 s | ~3 per combination | dt_crit; production runs use a fixed fraction (**OPEN**, default 0.7) |
| **D. Mesh-noise error bar** | 3-5 Halton seeds at the chosen setting | — | 3-5 | spread of the tip measures; a refinement level only counts as converged if its change exceeds this spread |
| **E. Actuation numerics** | `target_stiffness_factor` (e.g. 3, 7, 15, 30) × `STIFFNESS_RAMP_TIME` | short runs | ~6-8 | patch tracking error vs dt_crit; pick the factor |
| **F. Physics** | `FREQUENCY` and `AMPLITUDE` one at a time around the baseline (~5 levels each), then a ~4×4 grid over the interesting region; possibly `E_material` (0.45 vs 0.74 MPa, NOTES §4 issue 1) | at the settings from B-E | ~10 + ~16 | the study (**OPEN: ranges, levels, whether E is included**) |

### Stage A: tooling to build first

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

| # | decision | default suggestion |
|---|---|---|
| O1 | tolerance, and which tip measures (amplitude, phase lag, whole-trajectory RMS) | 2% on amplitude and phase |
| O2 | stiffness reference for stage B: Timoshenko on a single ribbon strip (as the wing's strip test), or a finite-element model of the whole fin | strip test first |
| O3 | `E_material` for the study: 0.74 MPa (no known source) or 0.45 MPa (UCL database, Dragon Skin 30) | 0.45 MPa |
| O4 | `mu` levels in stage C | 1, 3, 10 × water (30 × only if needed) |
| O5 | run length: cycles after the ramp | 3 |
| O6 | number of Halton seeds | 3 (5 if the spread is close to the tolerance) |
| O7 | fraction of dt_crit for production runs | 0.7 |
| O8 | baseline amplitude (2.5 or 3 cm) and the physics ranges/levels for stage F | 2.75 cm; f 0.25-1 Hz, A 1.5-3.5 cm, 5 levels |
| O9 | sweep-definition format for stage A (Python functions in `run_case.py`, or a config file per stage) | — |

## 7. Decisions log

Record each decision as it is made (date, decision, reason).

| date | # | decision | reason |
|---|---|---|---|
| 2026-10-18 | 1-5 | see §2 | planning discussion |
