# FinRay case: implementation notes

Context for anyone (or any new Claude session) picking up this case. Sections:

1. **Repository context**: where things live and how to run them
2. **Pipeline**: what gets called, in order
3. **Current model and why**: each design decision with its reason
4. **Open issues**: known problems and work not yet done, most important first
5. **History and verification log**: what changed when, superseded designs, and the checks behind the claims above

Last reorganised 2026-10-08. Sections 2-3 describe the code as it is; section 5 is history.

---

## 1. Repository context

- **Repo**: `fish-ib2d/` (local; not yet pushed). Moved here on 2026-10-08 from a fork of IB2d (`constancefr/IB2d`, `pyIB2d/Examples/FinRay`); all commits were carried over with rewritten paths.
- **Layout**: `cases/FinRay/` (this case), `cases/NACA0020_Wing/`, `common/mesh_utils.py` (shape-agnostic Halton fill + constrained Delaunay via meshpy/Triangle), `common/run_records.py` (run records + results index), `external/IB2d/` (upstream `nickabattista/IB2d`, git submodule pinned at `a5136ac`, **unmodified**).
- **Paths**: `main2d.py` and `FinRay_Geom.py` locate `external/IB2d/pyIB2d/IBM_Blackbox` and `common/` relative to their own file (`_REPO_ROOT` / `_COMMON_DIR`). Run every script **from inside `cases/FinRay/`**: IB2d reads `input2d` and `finray.*` from the current directory.
- **Python**: use `/opt/miniconda3/envs/ib2d/bin/python` (numpy, scipy, matplotlib, meshpy, numba, vtk, pyfftw). The system `python3` has no numpy. IB2d's solver requires numba.
- **Outputs**: `results/`, `viz_IB2d/`, `hier_IB2d_data/`, `*.mp4` are gitignored. The archived sweep results (~1.8 GB, 2026-09-30 to 10-01) were left in the old fork folder `ib2d_sim/IB2d/pyIB2d/Examples/FinRay/results/`.
- **Legacy**: `FinRay_Geom_zero_thickness.py` is the original single-point-thick implementation, kept for reference only.

## 2. Pipeline

**Sweep driver: `run_case.py`**
1. Loads `case_params.json` as `BASE_CASE`; `make_sweep()` varies `ribbon_fill_density_multiplier` over `[0, 1, 8, 32, 64]` × `ribbon_fill_stiffness_scale` over `[1.0]`.
2. For each case: writes `case_params.json`, patches `input2d` values by regex, runs `FinRay_Geom.py` then `main2d.py` **as subprocesses** (needed: `update_Target_Point_Positions.py` caches its baseline in module globals), moves outputs to a new folder `results/<case_name>_<YYYYmmdd-HHMMSS>/` (a rerun never overwrites an earlier run) with copies of `case_params.json` and `input2d`.
3. Writes that folder's **`run_info.json`**, the record of the run (`common/run_records.py`):
   - `inputs`: `geom` and `actuation` from `case_params.json`, plus every `input2d` value IB2d read.
   - `mesh`: from `mesh_info.json`, which `FinRay_Geom.py` writes each time it runs (gitignored): ds, fill spacing, multiplier, edge spacing, edge/fill points (fill requested and placed), points, triangles, springs, median/max k, k_Target, patch points, tail area, meshed-area fraction, build time.
   - `timing`: start/end, geometry and simulation wall time, steps, ms per step.
   - `outcome`: `completed` / `diverged` / `incomplete` / `no output`, first NaN time, last output time, exit code. A diverged case is flagged, not skipped.
   - `provenance`, captured at launch: commit, branch, IB2d submodule commit, package versions, machine. With uncommitted changes, the full diff (plus small untracked files) is saved as `code.diff` in the run folder, so quick tests before a commit stay reproducible.
4. After the sweep (even if a case crashed), rebuilds **`results/index.csv`**: one row per `run_info.json`, every field a column (`inputs.geom.E_material`, `mesh.n_fill_points`, `outcome.first_nan_time`, ...). Generated, never edited by hand: rebuild with `python ../../common/run_records.py index results`.
5. Runs made before records existed: `python ../../common/run_records.py backfill <run folder>...` writes a record from the folder's `case_params.json`/`input2d` and its output (mesh and code version unknown; start and wall time estimated from file modification times, dropped if the folder was copied). Done on 2026-10-13 for `results/`, `results_approach1/` and `old_results/FinRay_results_0/` (one folder there skipped: no `input2d`).

**Geometry: `FinRay_Geom.py` → `FinRay_Geom()`**
1. Parameters: grid (`Nx=32, Lx=1` → `dx=1/32`, `ds=dx/2=15.625 mm`), material and fill settings from `case_params.json["geom"]` (defaults in code).
2. `_Build_Tail_Mesh()` meshes the whole tail **once**:
   1. `_tail_outline()`: the tail as one polygon with 5 holes. `outer` = outer outline (a triangle); `windows` = the 5 open windows between base, rays and tip (4 trapezoids + 1 triangle). Built analytically from the half-width of the slanted edges moved by ±wall_thickness/2.
   2. Edge points: `mu.resample_loop()` puts points every ~h along `outer` and every window, keeping the corners, with `h = fill_spacing / sqrt(max(mult, 1))` = the fill's mean spacing.
   3. Fill: `mu.fill_polygon_with_halton(outer, inner_loops=windows)`, the whole tail in one call.
   4. `mu.triangulate_with_holes()`, once over all points (edge points first, then fill).
   5. `baseCornerIdx` = the 4 base corner points (outer and inner, bottom and top), found by coordinate.
   6. Returns `xLag, yLag, triangles, n_boundary, baseCornerIdx, outer, windows`.
3. `build_Tail_Connections()`: one spring per triangle edge, tributary-area `edge_k`; kind `'boundary'` (edge borders one triangle) or `'fill'`, for plotting only.
4. `_actuation_patch_indices()` → `k_Target` from `_max_incident_stiffness()`.
5. Preview plot → `geometry_preview.png` (edges, fill points, actuation patches as red stars).
6. Writes `finray.vertex`, `finray.geo_connect`, `finray.spring`, `finray.target`.

**Simulation: `main2d.py`** reads `input2d`, calls IB2d's `IBM_Driver`, which calls `update_Target_Point_Positions()` every step.

**Movies: `make_movie.py`** (VisIt CLI; `--all` renders every `results/` subfolder in one session).

**Fluid-free check: `verify_ribbon_fill.py`** builds the mesh at several densities and checks point-count scaling and degenerate springs. It does **not** calibrate stiffness.

### Parameters (`case_params.json`, current values)

| Key | Value | Meaning |
|---|---|---|
| `geom.E_material` | 0.74e6 Pa | Young's modulus (source unknown, see §4) |
| `geom.wall_thickness` | 3.96 mm | wall width (real build) |
| `geom.extrude_depth` | 0.027 m | real fin depth. **Not used by the simulation** (2D = per unit depth); only for scaling per-depth forces to 3D (F_3D ≈ F_2D · depth) |
| `geom.ribbon_fill_density_multiplier` | 64 | fill points = round(mult · tail area / ds²); edge points every ds/√mult |
| `geom.ribbon_fill_stiffness_scale` | 1.0 | multiplies **every** spring (calibration knob α) |
| `geom.ribbon_fill_spacing` | null → ds | Halton spacing |
| `geom.actuation_patch_radius` | null → wall_thickness | actuation patch radius |
| `geom.target_stiffness_factor` | 15 | k_Target / max incident spring stiffness |
| `actuation.FREQUENCY` | 0.1 Hz | (code default 1.0) |
| `actuation.AMPLITUDE` | 0.01 m | corner push in +x |
| `actuation.STIFFNESS_RAMP_TIME` | 0.2 s | target stiffness ramps 0 → full |
| `input2d` | `mu=1e-3, rho=1000 (water), dt=5e-4, Tfinal=15, Nx=Ny=32, Lx=Ly=1, supp=4` | springs=1, target_pts=1, update_target=1, beams=0 |

## 3. Current model and why

### 3.1 Geometry
- Triangle: length `L=0.136`, base width `W=0.051`, base at `x0=0.3`, centred on `y0=0.5`, tip pointing +x. Four rays parallel to the base at 25, 47, 69, 91 mm from it. All hard-coded in `_tail_outline` (TODO there: make these parameters).
- **The tail is one polygon with holes**: every wall (frame edges and rays) is `wall_thickness` wide, centred on its centerline, and the 5 open windows are the holes. Bending stiffness comes from the material in the wall (one side stretches while the other compresses) instead of a separately calibrated bending spring.
- **One mesh, no welding.** Rays and frame are meshed together, so the joints are meshed through like any other part of the wall. Before 2026-10-13 each ray was a separate ribbon welded to the frame by shared points (§5), which left each joint as a single ~4 mm edge with one triangle on each side.
- **Sharp corners**: the outline is the exact offset of the centerline triangle, so the outer tip is 10.7 mm past the centerline tip. The old ribbon code clamped the miter (≈7.9 mm). Tail area 1697.4 mm² (old: 1672.3).

### 3.2 Interior fill (Halton)
- Whole tail: `n = round(multiplier · tail_area / ds²)` fill points, i.e. multiplier 1 = one point per ds² of area. Edge points are spaced `h = ds/√multiplier` (ds at ×0/×1), the same as the fill's mean spacing, so edge and interior are equally dense:

  | multiplier | 0 | 1 | 4 | 8 | 16 | 32 | 64 | 128 |
  |---|---|---|---|---|---|---|---|---|
  | edge points | 51 | 51 | 105 | 146 | 211 | 301 | 420 | 599 |
  | fill points | 0 | 7 | 28 | 56 | 111 | 222 | 445 | 890 |
  | total Lagrangian points | 51 | 58 | 133 | 202 | 322 | 523 | 865 | 1489 |

- Mesh quality (×8 and up): no triangle under 5°, longest edge ≈ 1.8 h (3.5 mm at ×64). The old per-ribbon mesh had 7-8% of triangles under 5° and edges up to 23 mm, because its rails stayed at ~ds while the fill got denser.
- Fill points keep `min_boundary_gap` (≈0.5 mm) from every edge, so none lands almost on an edge segment.

### 3.3 Springs: one formula for all
Every spring is an edge of the tail's single constrained triangulation (edge + fill points; also at multiplier 0), with

```
k = to_ib(ribbon_fill_stiffness_scale · E · A_trib / L², ds)
  = ribbon_fill_stiffness_scale · E · A_trib / (L² · ds)
A_trib = 1/3 of the area of every triangle sharing the edge
```

- **Units** (`common/ib_units.py`): `E·A_trib/L²` is the physical stiffness per unit depth [N/m per m]. IB2d is 2D (everything per unit depth) and spreads `f·ds` with the global `ds = Lx/(2Nx)`, so the value written to `finray.spring` is divided by ds. `Nx`, `Lx` are read from `case_params.json["input2d"]`. Before 2026-10-08 the formula was `α·E·extrude_depth·A_trib/L²` (no ÷ds, ×0.05 m depth): 1280× too soft (§5).
- **Why**: this sets the springs' stored energy equal to the continuum strain energy of the triangles they span, so the wall's material is counted **exactly once at every density**. The previous scheme (rail springs, separate rung/diagonal cross springs, and fill springs on top) counted 175% (×1) to 234% (×64) of the real wall area, so the fin got stiffer as density rose and confounded any density study. Checked against the polygon areas (`verify_ribbon_fill.py`): exactly 100% at every density (×0-×256) since the switch to meshpy (§5).
- Outline springs border a triangle on one side only, so they get roughly half the k of a similar interior edge. The `'boundary'/'fill'` labels are for plotting only and don't change k.
- Spring rest length = the initial distance. `deg_NL = 1` (linear).
- Typical values (E = 0.74 MPa, scale 1, IB units): median k ≈ 1.0e7 from ×8 up (4e6-7e6 at ×0-×4); max 7e7-1.5e8 across densities.
- The triangulation is Triangle's constrained Delaunay (via `meshpy`, §5): it keeps every boundary segment, adds no points, and loses no wall area.

### 3.4 No beams
- `FinRay_Geom.py` writes no `.beam` file, and `input2d` has `beams=0` (true since the example's first commit).
- Why: bending comes from the triangulated ribbon. The old `.beam` file only held a weak (5% of k_Spring) per-rail regularisation against zig-zag, and IB2d never read it. A turning-angle check found no zig-zag (§5; that check was run on the older truss model). The NACA work also found that **IB2d's stock beam force is not conservative**: its left/right vertex signs are flipped, so it applies a net force. It can't be calibrated against EI, so it is best avoided altogether.

### 3.5 Actuation
- **Scheme** (`update_Target_Point_Positions.py`): two base-corner patches alternate. The top ("left") patch is pushed in +x by `AMPLITUDE · sin²(phase)` during the first half-cycle while the bottom patch holds still, then the roles swap. sin² gives zero velocity at both ends of each push. Target points are assigned to a patch by their rest y relative to the base midline.
- **Patch, not a single point**: pushing one point while its neighbour across the wall stayed rigid caused an unrealistic local twist. `_actuation_patch_indices` takes every **edge** point (index < `n_boundary`) within `actuation_patch_radius` of each corner, plus the 4 corner points. Fill points are excluded: forcing many of them into rigid lockstep caused NaN at high density.
- **Radius = wall_thickness, not ds**: a ds-sized radius reached a point on the *adjacent slant edge* (~13 mm away, closer than the base's own next point at ~17 mm) and caused NaN (×32 case). Since edge points now refine with density, the patch grows with it (same contact area): 4 points at ×0-×4, 8 at ×8-×16, 12-14 at ×32-×64, 20 at ×128.
- **Target stiffness**: `k_Target = target_stiffness_factor (15) × max over patch points of (sum of k of the springs attached to that point)`. With the single mesh: 6.8e8 at ×0/×1, 1.0e9-1.4e9 at ×4-×128 (IB units).
  - **Why relative, not absolute**: what matters is k_Target relative to the springs it must overpower. Too low and the patch lags the prescribed motion; too high and the explicit time stepping blows up. Both limits move with E, the stiffness scale, the density and ds.
  - **Why incident stiffness, not global max(edge_k)**: the global max swings ~3× between densities because of sliver triangles. The incident stiffness is steady: 4.5e7 at ×0/×1, 7.4e7-8.8e7 at ×4-×64 (IB units; ×1280 the pre-10-08 3.5e4 and 5.8e4-6.9e4).
  - **Why 15** (numbers in pre-10-08 units, all ×1280 now): it reproduces the previous effective value (1000·k_Spring × a runtime scale of 0.2 = 9.38e5 ≈ 13.5-16× incident at ×4-×64). The runtime `TARGET_STIFFNESS_SCALE` was removed so the stiffness is set in one place; the `.target` value is used as written, ramped in over `STIFFNESS_RAMP_TIME`.

## 4. Open issues (most important first)

1. **Stiffness calibration not built.** `scale = 1` is just the formula's literal value (now in correct units, §3.3). The planned check: a single ribbon, no fluid, one end clamped by stiff targets, a known tip force, compared with **Timoshenko** deflection `δ = FL³/(3EI) + FL/(κGA)`, κ ≈ 5/6, G = E/(2(1+ν)) with silicone's ν ≈ 0.499. Model it on `NACA0020_Wing/verify_structure.py`, which solves the static problem with IB2d's own force routines. Also: `E_material = 0.74 MPa` has no known source. The UCL Soft Robotics Materials Database (Dragon Skin 30, `Tensile-Tests-Data/Dragon Skin 30.csv`) gives E ≈ 0.44-0.45 MPa, near-linear to ~100% strain. No physical experiments are planned, so reference behaviour must come from analytical/FE models.
2. **Stability.**
   - **Long runs blow up at 8-14 s at every density** (15 s runs, dt 5e-4, old per-ribbon mesh). The fluid is the cause: its kinetic energy grows steadily from the start, almost all as grid-scale noise, centred just behind the actuated base corners, while spring strains stay at 2-7% until the last ~0.05 s. The cell Reynolds number (peak speed · dx / ν) is ~1000 after 1 s and ~6000 by 10 s; IB2d's central-difference advection went NaN at ~1250 on the NACA case. Fixes to try: raise `mu` (30× water was stable for NACA), refine the grid, or upwind advection (an IB2d change).
   - **×128 single mesh at dt 5e-4 goes NaN at t = 0.22 s**, just as k_Target finishes ramping in. Its 20-point actuation patch with k_Target per point is too stiff for this dt. Stable to 1 s with dt 2.5e-4, with the patch reduced to the 4 corner points, or at ×64 (14 points). Tip deflection at 1 s agrees within 0.3 mm across those three runs (`results/approach2_*_T1`). Open choice: lower dt at high density, cap the patch, or scale k_Target by 4/n_patch so the patch's total stiffness stays constant.
3. **Poisson's ratio ceiling.** A central-force triangular lattice has an effective ν ≈ 1/3 whatever the stiffness (Cauchy relation), while silicone has ν ≈ 0.5. The scale knob can only correct the effective Young's modulus.
5. **Triangles invert under large bending**, worst at the base. The inversions are small, localised and self-correcting; pure spring networks resist inversion but can't prevent it. Options if it matters: denser or stiffer fill at the base, or an anti-inversion penalty term.
6. **The patch radius is Euclidean, not edge-aware.** If the geometry changes (wall_thickness, rayPositions, L, W), re-check that the patch doesn't reach another edge, using the NaN check in §5.
7. **k_Target at ×0/×1** is about half the old value (5.2e5 vs 9.4e5, pre-10-08 units); check that the patch still follows the prescribed motion at low density.
8. Minor: the docstring of `update_Target_Point_Positions.py` says stiffness is "left untouched", but it is ramped. The geometry is hard-coded (4 rays; TODOs in `_Build_Tail_Geometry_Ribbon`). The TODO comment on `fill_density_multiplier` in `FinRay_Geom()` is answered in §3.2.

## 5. History and verification log

### Timeline
- **2026-09-09**: example created. Each edge was a single line of points with springs and beams, and rays were attached by an `'attach'` spring. That code is now in `FinRay_Geom_zero_thickness.py`.
- **2026-09-22**: two-rail ribbon model (rails + rung/diagonal "Warren truss" springs); Halton fill v0; actuation patch; automated movies.
- **2026-10-05**: single tributary-area spring formula (§3.3) replaces rail/rung/diagonal springs and `shear_Stiffness_From_Material`.
- **2026-10-06**: beams, `k_Spring`, `k_Beam` and the beam/stiffness helpers removed; k_Target rescaled (§3.5).
- **2026-10-08**: moved to `fish-ib2d`. The committed `finray.*` files were regenerated: they were left over from a ×32 sweep run and held the old k_Target (4.69e6, which would have been 5× too stiff once the runtime ×0.2 was removed). `run_case.py` regenerates them per case, so sweeps were never affected.
- **2026-10-08, units fix** (`common/ib_units.py`). IB2d is 2D, so every force is per unit depth, and it spreads `f·ds` with `ds = Lx/(2Nx)`. The old spring formula `α·E·extrude_depth·A_trib/L²` had no ÷ds (64× too soft) and multiplied by the 5 cm depth (20× too soft): together **1280× too soft**, and k_Target with it. Also `rho = 1` (water is 1000). Fixed: `k = α·E·A_trib/(L²·ds)`, rho = 1000, `Nx`/`Lx` read from `case_params.json` (they were hard-coded), and `extrude_depth` = 0.027 m (the real fin) kept for force reporting only. Every spring and k_Target is now exactly 1280× its old value. The wing case now uses the same helper, and its input files are byte-identical.
  - **Short runs** (2 s, scale 1, E 0.74 MPa): ×0 at dt 1e-3 is clean. ×64 at dt 1e-3 went NaN at t ≈ 0.22 s, just as k_Target finished ramping in; at dt 5e-4, 2.5e-4 and 1e-4 it is clean, with the same tip motion (Δ < 0.6 mm). **The default dt is now 5e-4.** Full-length stability is open issue 2.
  - **Tip motion barely changed**: at t = 2 s the tip moves −27.4 mm (×0) and −27.6 mm (×64), against −26.7 mm in the old-units ×0 run. At 0.1 Hz and 1 cm the motion is quasi-static, so fluid loads are negligible next to the elastic forces in both versions, and the shape is set by geometry plus the prescribed base motion. The units fix will matter for forces/thrust, at higher frequency or amplitude, and for stiffness calibration (issue 1).
- **2026-10-12, `mesh_utils` trimmed** (branch `approach1_trimmed`). Halton fill: scipy `qmc.Halton` replaces the hand-written sequence, a plain x/y sampling box replaces the PCA-oriented one, one distance helper replaces three, and the batch-doubling loop is gone. Same point counts and inside/gap rules, but different point positions. Triangulation: Triangle via `meshpy` replaces scipy Delaunay + the hand-written segment recovery (see the bug list), giving exactly 100% area at every density. `mesh_utils.py` went from 464 to 195 lines.
- **2026-10-13, single mesh** (branch `approach2`). The tail is now one polygon with 5 window holes (`_tail_outline`), with edge points every h = ds/√mult (`mu.resample_loop`), one Halton fill and one triangulation. This replaces the 5 separate ribbons, `_offset_polyline`, rail welding (`weld_to_inner_rail`, `mu.project_point_to_polyline`), `_n_points_for_length` and the per-ribbon spring merge. It closes the old open issue "thin edge triangles at high density". Effects: more points (×128: 1489 vs 954), the actuation patch grows with density, sharp tip (+1.5% area). `FinRay_Geom.py` went from 623 to 340 lines.
- **2026-10-13, run records**. Each run gets its own timestamped folder and a `run_info.json` (inputs, mesh metrics, timing, outcome, code version + `code.diff`), and `results/index.csv` is generated from those (§2). Older runs backfilled.

### Bugs found and fixed (all in git history)
- **Attach spring → weld**: removed the stiffest element; scale = 1.0 at ×1 went from NaN at t ≈ 3.3-3.6 s to clean through 5 s.
- **Rail desync**: inner-rail-only splicing broke the index pairing of the two rails, producing stray diagonal springs across the tail's open windows. Fixed by splicing a partner point into the outer rail.
- **Point-count rounding**: see §3.1. It also delayed the onset of triangle inversion (frame 0.44 → 0.62 s, ray0 0.9 → 1.58 s) without eliminating it.
- **matplotlib `Path(closed=True)` ignores the last vertex**, so every inside/outside test dropped one corner per loop (one ray lost ~70% of its area at ×16). Fixed by the `mu.closed_path()` helper, which appends the first vertex.
- **Unconstrained Delaunay dropped boundary segments** on the thin, curved frame (13-16 of 70 segments, ~9% of the area). `triangulate_with_holes` now recovers the missing segments itself (`_recover_loop_segments`, Anglada 1997): constrained Delaunay without a new dependency, since the `triangle` package has no wheel for Python 3.14. **Superseded 2026-10-12**: that recovery was buggy. It caused the old 0.2-0.5% area mismatch at high density and, with other fill positions, a dropped rail segment. One cause was found (the right-hand cavity chain was passed in reverse order), but fixing it still left 0.2-0.5% area errors, so there was at least one more bug, not traced. Replaced by Triangle via `meshpy`, which has a 3.14 wheel.
- **Lower-case keys**: `update_Target_Point_Positions.py` read lowercase `frequency` etc. while the JSON uses uppercase, so every actuation override silently did nothing.
- **Actuation drove only one rail point per corner**, because the target-row lookup (argmax/argmin) picked a single row; replaced by the patch.
- **Incomplete `holes → inner_loops` rename** caused a TypeError on every run, and a deleted `fill_spacing=None → ds` fallback broke `give_Me_Immersed_Boundary_Geometry`. Both fixed, and the output was checked to be unchanged (298 points at ×32).
- **Top-level sweep keys**: `run_case.py` keys must sit under `"geom"`; a top-level key was silently ignored.

### Superseded designs (for reading old commits)
- **Per-kind springs** (09-22 → 10-05): `k = EA_by_kind[kind] / L_own` per spring (rails E·A, rungs/diagonals G·A with ν = 0.499); fill springs used the tributary formula on top. That double-counted material and mixed two stiffness derivations at the rail boundary.
- **Beams**: `k_Beam = 0.05·k_Spring` per-rail regularisation, never active (beams=0). A real-EI version would need `EI/(ds·s⁵)` for IB2d's cross-product energy, but the stock beam force is non-conservative anyway (§3.4).
- **k_Spring** = `E·(t/2)·depth/ds`: from 10-05 to 10-06 it survived only to size k_Beam and k_Target.

### Checks and how to repeat them
- **NaN / divergence**: run to Tfinal, then look for `run_case.py`'s warning or `grep -i nan` the last `viz_IB2d/lagPtsConnect.*.vtk`. History: scale 1.0 and 0.5 at ×1 went NaN at t ≈ 3.3-3.6 s *before* the weld fix (0.1 and 0.01 were clean); clean after it. The ds-radius patch went NaN at ×32. The archived ×0-×64 runs at scale 1.0 are clean through 15 s (pre-10-05 model).
- **Zig-zag (beams needed?)**: signed turning angle at every rail point, rest vs t = 15 s (E 0.74 MPa, dt 1e-3, 0.1 Hz, 1 cm), at ×0 and ×32. Max change ~1.4-1.6°, no alternating-sign sawtooth. This was run on the **rung/diagonal truss model**, so re-run it if the current model shows artefacts or if AMPLITUDE or E change substantially.
- **Triangle inversion**: track each fill triangle's signed area against its rest sign over a run (×32, scale 1, 5 s). Frame peak 4/378 inverted, ray0 1/41 (after the rounding fix).
- **Material budget**: the sum over springs of `k·ds·L²/(E·scale)` should equal the total ribbon area (100%); `verify_ribbon_fill.py` checks it against the polygon areas.
- **Rigid patch**: at a given phase, every row in a patch gets the identical dx (`np.allclose`).
- **Mesh sanity**: `python verify_ribbon_fill.py` checks point-count ratio ≈ 4× per 4× multiplier and no degenerate springs.
