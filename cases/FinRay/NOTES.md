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
- **Layout**: `cases/FinRay/` (this case), `cases/NACA0020_Wing/`, `common/mesh_utils.py` (shape-agnostic Halton fill + constrained triangulation), `external/IB2d/` (upstream `nickabattista/IB2d`, git submodule pinned at `a5136ac`, **unmodified**).
- **Paths**: `main2d.py` and `FinRay_Geom.py` locate `external/IB2d/pyIB2d/IBM_Blackbox` and `common/` relative to their own file (`_REPO_ROOT` / `_COMMON_DIR`). Run every script **from inside `cases/FinRay/`**: IB2d reads `input2d` and `finray.*` from the current directory.
- **Python**: use `/opt/miniconda3/envs/ib2d/bin/python` (numpy, scipy, matplotlib, numba, vtk, pyfftw). The system `python3` has no numpy. IB2d's solver requires numba.
- **Outputs**: `results/`, `viz_IB2d/`, `hier_IB2d_data/`, `*.mp4` are gitignored. The archived sweep results (~1.8 GB, 2026-09-30 to 10-01) were left in the old fork folder `ib2d_sim/IB2d/pyIB2d/Examples/FinRay/results/`.
- **Legacy**: `FinRay_Geom_zero_thickness.py` is the original single-point-thick implementation, kept for reference only.

## 2. Pipeline

**Sweep driver: `run_case.py`**
1. Loads `case_params.json` as `BASE_CASE`; `make_sweep()` varies `ribbon_fill_density_multiplier` over `[0, 1, 8, 32, 64]` × `ribbon_fill_stiffness_scale` over `[1.0]`.
2. For each case: writes `case_params.json`, patches `input2d` values by regex, runs `FinRay_Geom.py` then `main2d.py` **as subprocesses** (needed: `update_Target_Point_Positions.py` caches its baseline in module globals), moves outputs to `results/<case_name>/` with copies of `case_params.json` and `input2d`.
3. Flags (does not skip) any case whose last `lagPtsConnect.*.vtk` contains NaN.

**Geometry: `FinRay_Geom.py` → `FinRay_Geom()`**
1. Parameters: grid (`Nx=32, Lx=1` → `dx=1/32`, `ds=dx/2=15.625 mm`), material and fill settings from `case_params.json["geom"]` (defaults in code).
2. `_Build_Tail_Geometry_Ribbon()`:
   1. Triangle centerline (base → top edge → tip → bottom edge), point counts from `_n_points_for_length()`.
   2. `_offset_polyline(closed=True)` → frame rails A/B; `mu.closed_path(...).contains_point` decides which is exterior.
   3. `baseCornerIdx` = 2 base corners × 2 rails (actuation seeds).
   4. For each of 4 rays: centerline → `_offset_polyline(closed=False)` → each of the 4 rail endpoints welded into the frame's inner rail by `weld_to_inner_rail()` (via `mu.project_point_to_polyline`; reuses a vertex within 2% of a segment end, otherwise splices a new point into the inner rail **and** a partner point into the outer rail) → ring polygon → `mu.fill_polygon_with_halton()` interior points.
   5. Frame interior: `mu.fill_polygon_with_halton(exterior, inner_loops=[inner])`.
   6. Returns `xLag, yLag, frameRailIdx, rayRailIdxList, baseCornerIdx, ribbonPolys`. Note `weld_to_inner_rail` mutates the frame rail index lists in place.
3. `build_Tail_Ribbon_Connections()`: rail-chain springs (dedupe self-loops/duplicates from welding) → `mu.triangulate_with_holes()` per ribbon → tributary area per edge → error if any rail segment isn't a triangulation edge → remaining edges added as `'fill'` → `edge_k`.
4. Preview plot → `geometry_preview.png`.
5. Writes `finray.vertex`, `finray.geo_connect`, `finray.spring`.
6. `_actuation_patch_indices()` → `k_Target` from `_max_incident_stiffness()` → `finray.target`.

**Simulation: `main2d.py`** reads `input2d`, calls IB2d's `IBM_Driver`, which calls `update_Target_Point_Positions()` every step.

**Movies: `make_movie.py`** (VisIt CLI; `--all` renders every `results/` subfolder in one session).

**Fluid-free check: `verify_ribbon_fill.py`** builds the mesh at several densities and checks point-count scaling and degenerate springs. It does **not** calibrate stiffness.

### Parameters (`case_params.json`, current values)

| Key | Value | Meaning |
|---|---|---|
| `geom.E_material` | 0.74e6 Pa | Young's modulus (source unknown, see §4) |
| `geom.wall_thickness` | 3.96 mm | ribbon width (real build) |
| `geom.extrude_depth` | 0.05 m | out-of-plane depth used in spring formula |
| `geom.ribbon_fill_density_multiplier` | 64 | interior points per ribbon = round(mult · area / ds²) |
| `geom.ribbon_fill_stiffness_scale` | 1.0 | multiplies **every** spring (calibration knob α) |
| `geom.ribbon_fill_spacing` | null → ds | Halton spacing |
| `geom.actuation_patch_radius` | null → wall_thickness | actuation patch radius |
| `geom.target_stiffness_factor` | 15 | k_Target / max incident spring stiffness |
| `actuation.FREQUENCY` | 0.1 Hz | (code default 1.0) |
| `actuation.AMPLITUDE` | 0.01 m | corner push in +x |
| `actuation.STIFFNESS_RAMP_TIME` | 0.2 s | target stiffness ramps 0 → full |
| `input2d` | `mu=1e-3, rho=1.0, dt=1e-3, Tfinal=15, Nx=Ny=32, Lx=Ly=1, supp=4` | springs=1, target_pts=1, update_target=1, beams=0 |

## 3. Current model and why

### 3.1 Geometry
- Triangle: length `L=0.136`, base width `W=0.051`, base at `x0=0.3`, centred on `y0=0.5`, tip pointing +x. Four rays parallel to the base at 25, 47, 69, 91 mm from it. All hard-coded in `_Build_Tail_Geometry_Ribbon` (TODO there: make these parameters).
- **Every edge is a 2-rail ribbon** of width `wall_thickness`, not a 1-point-thick line. This lets bending stiffness come from the material in the wall (one rail stretches while the other compresses) instead of a separately calibrated bending spring. `_offset_polyline` uses mitered normals, with the miter clamped (cos ≥ 0.25) so the tip doesn't spike.
- **Rays are welded to the frame** (shared vertices), not attached by a separate spring. The old `'attach'` spring (k = 50·k_Spring) was the stiffest element in the model and caused NaN. When a weld splices a point into the inner rail, it also splices a partner point into the outer rail; without it, springs that pair the two rails by array index connected the wrong points.
- **Boundary point counts** come from `_n_points_for_length`: round to the nearest number of *segments*, then points = segments + 1. Rounding the point count directly left base and ray segments at up to 214% of ds. Segments are now 70-124% of ds.

### 3.2 Interior fill (Halton)
- Per ribbon: `n = round(multiplier · ribbon_area / ds²)`. This answers the TODO in `FinRay_Geom()`: the multiplier scales the *area density*, and multiplier 1 means "one point per ds² of area". Because the ribbons are thin relative to ds, low multipliers give very few points:

  | multiplier | 0 | 1 | 4 | 8 | 16 | 32 | 64 |
  |---|---|---|---|---|---|---|---|
  | fill points | 0 | 6 | 27 | 56 | 109 | 220 | 439 |
  | total Lagrangian points | 78 | 84 | 105 | 134 | 187 | 298 | 517 |

  Ribbon areas: frame 1281.9 mm², rays 147.6 / 110.8 / 82.4 / 49.7 mm².
- The rail (boundary) spacing stays at ~ds regardless of the multiplier (see §4, thin edge triangles).

### 3.3 Springs: one formula for all
Every spring is an edge of one constrained triangulation per ribbon (rails + interior points; ribbons are triangulated even at multiplier 0), with

```
k = ribbon_fill_stiffness_scale · E · extrude_depth · A_trib / L²
A_trib = 1/3 of the area of every triangle sharing the edge (summed across ribbons)
```

- **Why**: this sets the springs' stored energy equal to the continuum strain energy of the triangles they span, so the wall's material is counted **exactly once at every density**. The previous scheme (rail springs, separate rung/diagonal cross springs, and fill springs on top) counted 175% (×1) to 234% (×64) of the real wall area, so the fin got stiffer as density rose and confounded any density study. Checked: 100% at ×0, 1, 4, 16, 64, 256.
- Rail springs border triangles on one side only, so they get roughly half the k of a similar interior edge. The `'frame'/'ray'/'fill'` labels are for plotting only and don't change k.
- Spring rest length = the initial distance. `deg_NL = 1` (linear).
- Typical values (E = 0.74 MPa, scale 1): median k ≈ 5-7e3, max 5e4-1.5e5 depending on density. The max is set by whichever thin "sliver" triangle the fill happens to create, so it is noisy.
- Two triangulation bugs were fixed in `mesh_utils` (§5). The triangulation now loses no wall area and drops no boundary segments.

### 3.4 No beams
- `FinRay_Geom.py` writes no `.beam` file, and `input2d` has `beams=0` (true since the example's first commit).
- Why: bending comes from the triangulated ribbon. The old `.beam` file only held a weak (5% of k_Spring) per-rail regularisation against zig-zag, and IB2d never read it. A turning-angle check found no zig-zag (§5; that check was run on the older truss model). The NACA work also found that **IB2d's stock beam force is not conservative**: its left/right vertex signs are flipped, so it applies a net force. It can't be calibrated against EI, so it is best avoided altogether.

### 3.5 Actuation
- **Scheme** (`update_Target_Point_Positions.py`): two base-corner patches alternate. The top ("left") patch is pushed in +x by `AMPLITUDE · sin²(phase)` during the first half-cycle while the bottom patch holds still, then the roles swap. sin² gives zero velocity at both ends of each push. Target points are assigned to a patch by their rest y relative to the base midline.
- **Patch, not a single point**: pushing one point while its neighbour across the wall stayed rigid caused an unrealistic local twist. `_actuation_patch_indices` takes every **boundary** point within `actuation_patch_radius` of each corner. Interior fill points are excluded: forcing many of them into rigid lockstep caused NaN at high density.
- **Radius = wall_thickness, not ds**: a ds-sized radius reached a point on the *adjacent slant edge* (~13 mm away, closer than the base's own next point at ~17 mm) and caused NaN (×32 case). At the default radius the patch is exactly the 4 corner rail points at every density.
- **Target stiffness**: `k_Target = target_stiffness_factor (15) × max over patch points of (sum of k of the springs attached to that point)`. At ×64 this gives 1.04e6.
  - **Why relative, not absolute**: what matters is k_Target relative to the springs it must overpower. Too low and the patch lags the prescribed motion; too high and the explicit time stepping blows up. Both limits move with E, the stiffness scale, the density and ds.
  - **Why incident stiffness, not global max(edge_k)**: the global max swings ~3× between densities because of sliver triangles. The incident stiffness is steady: 3.5e4 at ×0/×1, 5.8e4-6.9e4 at ×4-×64.
  - **Why 15**: it reproduces the previous effective value (1000·k_Spring × a runtime scale of 0.2 = 9.38e5 ≈ 13.5-16× incident at ×4-×64). The runtime `TARGET_STIFFNESS_SCALE` was removed so the stiffness is set in one place; the `.target` value is used as written, ramped in over `STIFFNESS_RAMP_TIME`.

## 4. Open issues (most important first)

1. **Units and force scaling: audit before trusting any magnitudes.**
   - *Verified*: IB2d treats every Lagrangian force as a density and spreads `f·ds` with the global `ds = Lx/(2Nx)` (`please_Find_Lagrangian_Forces_On_Eulerian_grid.py:316`). A spring constant k therefore acts physically as `k·ds`, i.e. **~64× softer** than the formula's value at Nx=32. The NACA case corrects for this (`K_phys/ds`); FinRay does not. Ratios such as k_Target / incident stiffness are unaffected.
   - *Not verified, worth checking in the same audit*: IB2d is 2D, so its forces are per unit depth, yet the spring formula multiplies by `extrude_depth` (0.05 m). Also, `input2d` has `rho = 1.0` with `mu = 1e-3`; water would be rho = 1000.
2. **Stiffness calibration not built.** `scale = 1` is just the formula's literal value. The planned check: a single ribbon, no fluid, one end clamped by stiff targets, a known tip force, compared with **Timoshenko** deflection `δ = FL³/(3EI) + FL/(κGA)`, κ ≈ 5/6, G = E/(2(1+ν)) with silicone's ν ≈ 0.499. Model it on `NACA0020_Wing/verify_structure.py`, which solves the static problem with IB2d's own force routines. Also: `E_material = 0.74 MPa` has no known source. The UCL Soft Robotics Materials Database (Dragon Skin 30, `Tensile-Tests-Data/Dragon Skin 30.csv`) gives E ≈ 0.44-0.45 MPa, near-linear to ~100% strain. No physical experiments are planned, so reference behaviour must come from analytical/FE models.
3. **The current model has never been run to full length.** Every archived result in `results/` (all NaN-free through 15 s at scale 1.0, ×0-×64) dates from 2026-09-30 to 10-01, **before** the single-formula springs (10-05) and the k_Target change (10-06). The stability and zig-zag findings in §5 were made on older models too. Re-run the sweep, and repeat the turning-angle zig-zag check if anything looks off.
4. **Poisson's ratio ceiling.** A central-force triangular lattice has an effective ν ≈ 1/3 whatever the stiffness (Cauchy relation), while silicone has ν ≈ 0.5. The scale knob can only correct the effective Young's modulus.
5. **Thin edge triangles at high density.** Rail points stay ~ds (15 mm) apart while interior spacing shrinks (~2 mm at ×64), giving long, thin boundary triangles and sliver springs. Consider refining the rails along with the fill density.
6. **Triangles invert under large bending**, worst at the base. The inversions are small, localised and self-correcting; pure spring networks resist inversion but can't prevent it. Options if it matters: denser or stiffer fill at the base, or an anti-inversion penalty term.
7. **The patch radius is Euclidean, not edge-aware.** If the geometry changes (wall_thickness, rayPositions, L, W), re-check that the patch doesn't reach another edge, using the NaN check in §5.
8. **k_Target at ×0/×1** is about half the old value (5.2e5 vs 9.4e5); check that the patch still follows the prescribed motion at low density.
9. Minor: the docstring of `update_Target_Point_Positions.py` says stiffness is "left untouched", but it is ramped. The geometry is hard-coded (4 rays; TODOs in `_Build_Tail_Geometry_Ribbon`). The TODO comment on `fill_density_multiplier` in `FinRay_Geom()` is answered in §3.2.

## 5. History and verification log

### Timeline
- **2026-09-09**: example created. Each edge was a single line of points with springs and beams, and rays were attached by an `'attach'` spring. That code is now in `FinRay_Geom_zero_thickness.py`.
- **2026-09-22**: two-rail ribbon model (rails + rung/diagonal "Warren truss" springs); Halton fill v0; actuation patch; automated movies.
- **2026-10-05**: single tributary-area spring formula (§3.3) replaces rail/rung/diagonal springs and `shear_Stiffness_From_Material`.
- **2026-10-06**: beams, `k_Spring`, `k_Beam` and the beam/stiffness helpers removed; k_Target rescaled (§3.5).
- **2026-10-08**: moved to `fish-ib2d`. The committed `finray.*` files were regenerated: they were left over from a ×32 sweep run and held the old k_Target (4.69e6, which would have been 5× too stiff once the runtime ×0.2 was removed). `run_case.py` regenerates them per case, so sweeps were never affected.

### Bugs found and fixed (all in git history)
- **Attach spring → weld**: removed the stiffest element; scale = 1.0 at ×1 went from NaN at t ≈ 3.3-3.6 s to clean through 5 s.
- **Rail desync**: inner-rail-only splicing broke the index pairing of the two rails, producing stray diagonal springs across the tail's open windows. Fixed by splicing a partner point into the outer rail.
- **Point-count rounding**: see §3.1. It also delayed the onset of triangle inversion (frame 0.44 → 0.62 s, ray0 0.9 → 1.58 s) without eliminating it.
- **matplotlib `Path(closed=True)` ignores the last vertex**, so every inside/outside test dropped one corner per loop (one ray lost ~70% of its area at ×16). Fixed by the `mu.closed_path()` helper, which appends the first vertex.
- **Unconstrained Delaunay dropped boundary segments** on the thin, curved frame (13-16 of 70 segments, ~9% of the area). `triangulate_with_holes` now recovers the missing segments itself (`_recover_loop_segments`, Anglada 1997): constrained Delaunay without a new dependency, since the `triangle` package has no wheel for Python 3.14.
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
- **Material budget**: the sum over springs of `k·L²/(E·depth·scale)` should equal the total ribbon area (100%).
- **Rigid patch**: at a given phase, every row in a patch gets the identical dx (`np.allclose`).
- **Mesh sanity**: `python verify_ribbon_fill.py` checks point-count ratio ≈ 4× per 4× multiplier and no degenerate springs.
