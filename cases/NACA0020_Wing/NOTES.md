# NACA0020 wing case: implementation notes

Replicates the soft NACA0020 wing of `WaterLily_FSI/WaterLily_2Way_Coupling_Example/gvs_wings/Wing_3order.jl`, as run by `gvs_tests.jl` (the "SoRoLily" model, validated against experiments: **the ground truth**). Built with the same method as `cases/FinRay` ("approach 2": outline the shape, fill it, triangulate once, one spring per triangle edge). Sections:

1. **Repository context**: where things live and how to run them
2. **Pipeline**: what gets called, in order
3. **Current model and why**: the reference model, the mapping, the calibration
4. **Open issues and guesstimates**: known problems, most important first
5. **History**: the previous (centreline) model and its validation runs

Last reorganised 2026-10-15; §1-2 updated 2026-10-18 (cluster tooling, hooks, reference extraction).

---

## 1. Repository context

- Same layout and conventions as `cases/FinRay/NOTES.md` §1: run every script **from inside `cases/NACA0020_Wing/`**, with the `ib2d` conda env. IB2d is the unmodified `external/IB2d` submodule. On the cluster: `slurm/submit_sweep.sh` (as FinRay, `CLUSTER.md` §4; extra `TIME_LIMIT`, `MEM`), results in `~/finray_fsi/sim_data/fish-ib2d/NACA0020_Wing/<sweep>/`.
- **`IMPLEMENTATION_DECISIONS.md`**: every modelling decision against SoRoLily, verified at source, with the defaults adopted; **`RESULTS.md`**: what the runs showed.
- Shared code: `common/mesh_utils.py`, `common/ib_units.py`, `common/run_records.py` (`run_case()`, run records, index).
- `wing_case.py` reads `case_params.json` and `input2d` (the wing reads the grid from `input2d`; `run_case` writes the same values to both).
- Reference: the SoRoLily source is in `~/finray_fsi/WaterLily_FSI/WaterLily_2Way_Coupling_Example/`. `reference/extract_reference.py` turns its `data_gvs/gvs_tests.h5` into the committed `reference/sorolily_reference.npz` (physical units; TE, centreline, Fx, Fy rebuilt from the panel pressures, surface pressures); `reference/sorolily.py` ports its kinematics, panels and filters. `force_history[2]` (Julia indexing) is the summed **x**-force although it is plotted as "Lift"; total Fy is not stored.
- Legacy: `Wing_Geom_centerline.py` (+ the legacy branch of the user-force hook) = the previous model (§5); `analyze_run.py` = the previous post-processing, superseded by `postprocess.py`.

## 2. Pipeline

```sh
python Wing_Geom.py            # wing.* files, wing_mesh.npz, mesh_info.json, geometry_preview.png
                               #   (stiffness_scale: null -> measured for this mesh, ~20-60 s)
python verify_structure.py     # the full fluid-free calibration study (density sweep), verify_structure.png
python spring_force_check.py   # the hook's vectorised springs == IB2d's (after Wing_Geom.py)
python run_case.py --sweep NAME [--list | --case N]   # geometry + FSI run + postprocess.py
python report.py SWEEP_FOLDER  # report/summary.md + figures for a cluster sweep
```

- **`Wing_Geom.py`**: `_wing_outline` → `build_wing_mesh` (Halton fill + one constrained triangulation) → `mu.tributary_springs` → targets → `stiffness_scale` measured if null (`verify_structure.measure_scale`) → `bend` (pre-bent initial state) → files. Also writes `wing.user_force` (the hook's mode) and `wing_mesh.npz` (rest state, triangulation, clamp: read by the hooks and `postprocess.py`).
- **Hooks** (IB2d unchanged): `give_Me_General_User_Defined_Force_Densities.py` (`user_force_model = 1`) computes the springs vectorised (`geom.spring_force = "hook"`, `input2d springs = 0`) and records every Lagrangian point and the clamp force (`viz_IB2d/wing_history.npz`: points every 0.01 s, force every step); `please_Compute_External_Forcing.py` is the free stream (§3.4) + start-up ramp, and records the fluid's mean velocity, KE and grid-scale KE (`viz_IB2d/fluid_history.npz`). `wing_diagnostics.py` is their shared recorder.
- **`run_case.py`** / **`sweeps.py`**: as FinRay (`run_records.run_case()` per case; sweeps are functions in `sweeps.py`). Each run folder gets `run_info.json`, then **`postprocess.py`** adds the comparison with SoRoLily (`analysis`), `comparison.png`, `flow_snapshots.png`, `timeseries.csv`. Crashed runs are post-processed too (IB2d stops with an IndexError on NaN or overflowed positions rather than writing NaN frames; `analysis.stable` is the reliable flag).

### Parameters (`case_params.json`)

| Key | Value | Meaning |
|---|---|---|
| `geom.E_flex` | 1e4 Pa | Young's modulus of the flexible link (reference) |
| `geom.rigid_length` | 0.06 m | clamped leading-edge block (reference link 2) |
| `geom.init_q` | [1, 0, 0, 0] | initial curvature modes (reference `init_state`) |
| `geom.fill_density_multiplier` | 4 | fill points = round(mult · area / fill_spacing²); edge points every h = fill_spacing/√mult |
| `geom.fill_spacing` | null → ds | Halton spacing (so h/dx is the same at every `Nx`) |
| `geom.stiffness_scale` | **null** | measured for each mesh at geometry time (4.083 at ×4, Nx 192; §3.3) |
| `geom.target_stiffness_factor` | 3 | k_Target / max incident spring stiffness; 15 halves dt_crit (`IMPLEMENTATION_DECISIONS.md` O4) |
| `geom.spring_force` | "hook" | springs computed in the user-force hook (vectorised) |
| `flow.*` | U = 0.2 m/s, AoA 0° → 15° at 5 s over 0.1 s, fringe 40 s⁻¹, mean 20 s⁻¹, `start_ramp` 0.02 s | free stream (§3.4) |
| `output.field_interval` | 0.05 s | `print_dump` = this / dt (vorticity + pressure only) |
| `diagnostics.interval`, `save_every` | 0.01 s, 1 s | hook sampling, and how often the histories are written |
| `input2d` | rho 1000, dt 2e-4, Tfinal 10, 192², 1.2 m box, mu 30× water | production dt per Nx: `sweeps.DT_PRODUCTION` |

## 3. Current model and why

### 3.1 What the reference model is

Superseded by `IMPLEMENTATION_DECISIONS.md` §3, which re-checked every item at source (2026-10-18) and corrects several (fixed 0.01 s coupling step, filtered staggered coupling, boundary and initial conditions, QUICK advection). Kept as written on 2026-10-15:

| GVS (`Wing_3order.jl` + `gvs_functions.jl`) | Consequence |
|---|---|
| link 1 = massless carriage, `P`+`R` joints, but `dynamics_gvs` sets `state_dot[7:8] = 0` and `init_state` has zero velocity | heave and pitch are **frozen**: the wing is a fixed body |
| link 2, 0.06 m, `Phi_dof = 0` → rigid | fixed leading-edge block, 0–0.06 m from the LE (30% chord) |
| link 3, 0.14 m, curvature only, cubic Legendre basis (`Phi_dof=(1,0,0)`, `Phi_odr=(3,0,0)`) | **inextensible, unshearable** cantilever clamped at 0.06 m, whose curvature is a cubic polynomial (4 numbers) |
| `E = 1e4 Pa`, `ρ = 1064.6`, `η = 2e3 Pa s`, `ν = 0.5` | EI(x) = E·h·w(x)³/12; Kelvin-Voigt bending damping (not reproduced, §4) |
| `w(x)` fitted thickness, zero at the tip, ≈ 0.04 m at x = 0.06 m | NACA0020-like, chord 0.2 m. The fit crosses zero 0.09 mm before the TE (it is −4e-5 m at x = 0.2; `thickness()` clips it) |
| `h = 0.1 m` span, forces × `w = 0.1` | structure and fluid load both ∝ span ⇒ **IB2d uses span = 1 (per unit depth)** |
| tendons `n_sact = 2` | `dynamics_gvs` multiplies the tension by 0 → **no actuation** |
| `init_state = [0,0,1,0,0,0, 0…]` | starts **pre-bent**: κ = 1/L₃ = 7.14 m⁻¹, tip rotated 1 rad, TE 64.4 mm up, released at t = 0 |
| inflow U = 0.2 m/s, α: 0° → 15° at t = 5 s (0.1 s smoothstep), ν = 1.0023e-6, ρ = 1000; pivot (0.4, 0.6) m in a 1.2 × 1.2 m box; 32 cells / chord | reproduced (fringe forcing, 192², same pivot) |
| the fluid sees a degree-1 NURBS (a polygon) through ±w/2 at ~16 stations per side (`body_simulations.jl:134-140`) | IB2d uses the smooth fit instead; the difference is under one grid cell |

### 3.2 Geometry and mesh (`Wing_Geom.py`)
- **One polygon, one mesh.** `_wing_outline(h)`: LE → upper surface → TE → lower surface, points every ~h in arc length (`surface_stations`, refined near the nearly vertical LE), with a vertex on each surface exactly at x′ = 0.06. Halton fill (`mu.fill_polygon_with_halton`), one Triangle triangulation (`mu.triangulate_with_holes`). Mesh area = polygon area (100.0000%).
- **Clamp line.** The section x′ = 0.06 m is a chain of points (spacing h) passed to Triangle as extra edges (`extra_segments`), so no triangle crosses it and the flexible part is clamped exactly at 0.06 m, as in the reference. Fill points keep their gap from it (`keep_out_segments`).
- **`min_boundary_gap = h/2`**, passed explicitly. `mu`'s default (0.25·area/perimeter) would be ~3.3 mm here, more than ds and than the TE's thickness, which would leave a fill-free band along the skin.
- **Densities** (ds = 3.125 mm at 192², 1.2 m):

  | multiplier | 1 | 4 | 8 | 16 | 32 |
  |---|---|---|---|---|---|
  | points | 703 | 2520 | 4870 | 9498 | 18655 |
  | edge spacing h [mm] | 3.1 | 1.6 | 1.1 | 0.78 | 0.55 |

  At ×8: 413 edge points, 4457 fill, min triangle angle 6.2°, no triangle under 5°.

### 3.3 Springs and stiffness calibration
- Every triangle edge is a spring, `k = stiffness_scale · E · A_trib / (L² · ds)` (`mu.tributary_springs`, the FinRay formula), rest length = the **straight** wing's edge length.
- **Why `stiffness_scale` is ~4, not 1.** The tributary formula puts each triangle's area into its springs exactly once, but that does not make the network as stiff as a sheet of modulus E. A spring only resists stretching along its own direction, so when a triangle is stretched one way its three edges store only part of the energy a real solid would. Measured on Halton + Triangle meshes (uniform strain, 2026-10-15): effective E = 0.33 E (×1), 0.36 E (×8, ×64); exactly E/3 for an equilateral mesh. Effective Poisson ratio 0.29–0.34 (§4). **This applies to FinRay too** (its open issue 1): `scale = 1` there is ~3× too soft.
- **`verify_structure.py`** (fluid-free, linear static solve with the sparse stiffness matrix of the springs and targets):
  1. *strip*: a 0.14 × 0.007 m strip (L/t = 20, shear < 0.3%), clamped, uniform load, vs Euler–Bernoulli → the pure network correction: 2.42 (×1), 3.14 (×4), 3.35 (×8), 3.62 (×16), 3.65 (×32). It rises with density because a finer random network has more ways to deform locally.
  2. *wing*: uniform load 0.5 N/m on the flexible upper skin; the centreline deflection (interpolated on the rest triangulation) is fitted by least squares to **the reference's own model** (Ritz with its cubic curvature basis) over 0 < s/L₃ ≤ 0.9 → 3.74 (×1), 4.07 (×4), **3.97 (×8)**, 4.05 (×16), 3.98 (×32). Fitting to Euler–Bernoulli with the true thickness instead gives the same within 0.1%, because the two agree to 0.2% up to s/L₃ = 0.9 (table below).
  3. *IB2d cross-check*: IB2d's own spring + target routines at the solved shape balance the load. The residual / load ratio falls exactly with the load (13, 0.13, 0.0013 at load ×1, ×1e-2, ×1e-4), so the linear part matches exactly. What remains is the rotation term a linear solve leaves out, which is large next to the load here only because the wing is far stiffer in stretching than in bending.
  4. *clamp*: target_stiffness_factor 150 instead of 15 changes the tip deflection by 0.04%.
- **Why the wing needs more than the strip** (3.97 vs ~3.4–3.6): shear. The reference cannot shear, but the meshed solid does, and near the root the link is short and thick (L₃/w ≈ 3.5). The fitted profile shows it: the mesh deflects more than the reference near the clamp (shear grows ∝ s there, bending ∝ s²) and slightly less near the tip (table). Fitting to the reference puts that into `stiffness_scale` (9% above the strip's value at ×32, 18% at ×8), as intended: the goal is to match the reference's behaviour, not a shearable solid's.
- **Why the fit stops at s/L₃ = 0.9.** The last 10% is where the TE is under ~9 mm thick and the reference's cubic curvature cannot follow the real (or meshed) thin edge. Up to 0.9 the reference agrees with beam theory with the true thickness to 0.2%, so this is the reference model's limit, not a material difference:

  | s/L₃ | 0.12 | 0.25 | 0.5 | 0.75 | 0.87 | 1.0 (tip) |
  |---|---|---|---|---|---|---|
  | mesh ×8, scale 3.966 [mm] | 0.0173 | 0.0582 | 0.1943 | 0.4120 | 0.5633 | 0.7724 |
  | reference (cubic curvature) | 0.0131 | 0.0501 | 0.1866 | 0.4128 | 0.5744 | 0.7855 |
  | beam, true thickness | 0.0132 | 0.0500 | 0.1867 | 0.4128 | 0.5727 | 0.8336 |

  Tip / reference tip: 0.94 (×1), 0.98 (×4), 0.98 (×8), 0.99 (×16, ×32).
- **Default density ×8**: its scale is within 0.5% of ×32's. The ×4 to ×32 values scatter by ±1.5% with no trend (each density is a different Halton sample), so **`stiffness_scale` belongs to one mesh**: re-run `verify_structure.py` and update `case_params.json` whenever the multiplier, fill_spacing or grid changes.
- The old NOTES' "EB tip 4.9 mm, log-divergent without a floor" was an artefact of the fit's zero-thickness sliver at the TE (above) with a 1 µm floor on a coarse grid. With a 10 µm floor it converges to 0.834 mm. The triangulated wing uses the true shape, so `EI_floor_thickness` is gone.

### 3.4 Clamp, initial state, free stream
- **Clamp**: every point with x′ ≤ 0.06 m (rigid block + clamp line) is a target point, `k_Target = target_stiffness_factor × max incident spring stiffness` (FinRay rule; factor 15 until 2026-10-15, **3 since 2026-10-18**: 15 halves dt_crit and changes the static tip by 0.04 %). `update_target = 0`: the block never moves. Springs with both ends clamped never stretch, so they are left out of `wing.spring` (5068 of 14229 at ×8). They are still in `wing.geo_connect` and `geometry_preview.png`, but **VisIt does not show them**: IB2d draws `lagPtsConnect.*.vtk` from the springs (`IBM_Driver.py:658`, `give_Me_Lag_Pt_Connects`, lines 1321-1346) and never reads `.geo_connect`, so the clamped block appears as points without edges. **This is needed for stability, not just speed**: with them, the max incident stiffness (and so k_Target) was 4.5× higher (6.3e8 vs 1.4e8 IB units at ×8) and the run went NaN at t = 0.002 s; without them it is stable (§5, 2026-10-15).
- **Pre-bend** (`bend`): every flexible point (x′, y′) → c(s) + y′·n(s), s = x′ − 0.06, on the bent centreline of `init_q` (plane sections stay plane and normal, as in the reference's beam kinematics). `.vertex` holds the bent shape, `.spring` the straight rest lengths, so the bend is the initial elastic state. TE starts 64.36 mm up (reference 64.4 mm).
- **Start-up** (2026-10-18): the mean flow is driven 0 → U along a smoothstep in `flow.start_ramp` = 0.02 s, standing in for WaterLily's u = U at t = 0. A one-step impulse instead excites the clamp's ringing mode (the wing's interior fluid on the stiff tethers, period ~1.8 ms) to ±4000 N/m. Some ringing (~20-30 N/m RMS at ~600 Hz) is always present in the tether force; it is filtered out of the force comparison and moves the fluid by ~1 mm/s.
- **Free stream** (`please_Compute_External_Forcing.py`): IB2d is doubly periodic, so a fringe at x ∈ [1.0, 1.2] relaxes u → U(cos α, sin α) at 40 s⁻¹, plus a domain-mean-velocity controller (without it the turned flow reached the wing ~2 s after the AoA switch; with it, within ~0.2 s).

## 4. Open issues and guesstimates (most important first)

1. **Viscosity floor (2026-10-18, `RESULTS.md` §3).** The schedule completes only at ≥ 30× water (Re 1330), at Nx 192 and 384; everything below diverges in 0.2-2 s through a grid-scale (odd-even) checkerboard of IB2d's collocated central-difference solver. Not a dt limit, not the release, not the structural mesh. Water would need ~1000 cells/chord, if it is reachable at all.
2. **Bistable deflected wing (2026-10-18, `RESULTS.md` §4).** At 30× the wing settles at ±22-29 mm TE deflection at 0° (side set by the pre-bend's sign) and stays deflected at 15°; the reference: 0.3 / 4.9 mm. Robust to fringe, start-up, clamp, mesh, grid and Re 400-1330; the flap carries +0.5 N/m against −3.6 N/m on the block (static divergence). Separating the Re floor from the missing damping needs a SoRoLily re-run (ν × 30, η = 0). Even at 30× the straight start and a half pre-bend diverge.
3. **No Kelvin-Voigt damping** (unchanged from the old model). The reference is overdamped (3η/E ≈ 0.6 s): its wing filters the flow loads. An explicit dashpot on IB2d's massless points has loop gain c·dt·g/m ~ 1e2–1e3 (nodal inertia ~ρ·dx²), so it is unstable or limited to ~1% of η. Needs implicit coupling or a partitioned scheme; this is the main expected dynamic difference.
4. **Run time / dt (resolved 2026-10-18).** Vectorised springs in the hook: 24 ms/step at Nx 192 (20 min per 10 s run), ~150 ms at Nx 384 (4 h), cluster batch nodes. dt_crit (factor 3) = 1e-3, 3.125e-4, 1.5625e-4 s at Nx 96, 192, 384.
5. **Trailing edge**: thickness → 0, so the last few mm are one or two triangles thick. The meshed TE can curl; the reference's cannot (tip +6% for true-thickness beam theory under a uniform load). Expect TE-deflection comparisons to show it.
6. **Poisson's ratio 0.29–0.34, not 0.5**: a central-spring network can't do better (FinRay issue 3). Little effect on the bending of a slender link; some on how the thick root section deforms.
7. **Guesstimates** (2026-10-18: factor now 3; fringe ×0.5/×2 changes the TE means by ≤ 10 %, `stage5_sensitivity`). Before: `target_stiffness_factor` 15 is FinRay's value (insensitive here, §3.3); fringe rate 40 s⁻¹ and mean rate 20 s⁻¹ were set on the old model (only the AoA arrival time was checked); the calibration load (uniform per length on the upper skin) is not the flow's pressure distribution, but the fit is linear and dominated by the bending of the root half.
8. **Density**: the airfoil interior is fluid at 1000 kg/m³ that moves with the wing, vs ρ_s = 1064.6 (+6%). No mass points.
9. **Fill inside the clamped block** (~1500 target points at ×8) does nothing mechanically. Leaving it out (fill only the flexible part) would save points but makes the mesh non-uniform; not done.

## 5. History

### Timeline
- **2026-09** (commit 280cce2): case created as the **centreline model** (`Wing_Geom_centerline.py`): a 10 mm centreline chain with a turning-angle bending energy ½(EI/s)θ² per vertex through `user_force_model` (IB2d's stock `.beam` is not the gradient of its energy, so it can't be calibrated to EI; see `give_Me_General_User_Defined_Force_Densities.py`), axial springs for inextensibility, a skin of outline points tied to the chain by two stiff springs each, an 8 mm thickness floor on EI, and absolute `K_target = 1e4`. Calibration: tip 0.772 mm vs EB 0.770 (+0.2%, load on the centreline), 0.74 (−4%, load on the skin).
- **2026-10-08**: moved to `fish-ib2d`; units helper `common/ib_units.py` (input files byte-identical).
- **2026-10-15**: rebuilt with the FinRay method (§3). Shared code moved to `common/` (`tributary_springs`, `max_incident_stiffness`, `run_case`, `patch_input2d`); FinRay's generated files are byte-identical at ×1 and ×64. `verify_structure.py` rewritten for the triangulated wing; `analyze_run.py` reads the run folder's `wing.target` and `run_info.json`.
  - First FSI smoke run (0.5 s, ×8, dt 1e-4) went NaN before t = 0.05 s. Short scratch runs (0.02-0.06 s) were all clean once the springs between clamped points were dropped; keeping them reproduced the NaN at t = 0.002 s (k_Target 6.3e8 vs 1.4e8, §3.4). The repeated 0.5 s run (`results/smoke2_0p5s_*`) completed: TE 64.4 → 54 → 31 → 19 → 11 → 4 mm at t = 0, 0.1, 0.2, 0.3, 0.4, 0.5 s (reference: 30 mm at 0.5 s, slowed by its KV damping; the old model was also much faster). Release-transient forces of 10-30 N/m (sampled every 0.05 s).
- **2026-10-18**: SoRoLily re-read at source (`IMPLEMENTATION_DECISIONS.md` §3); cluster tooling (`sweeps.py`, `slurm/`, `postprocess.py`, `report.py`, `summarise.py`), hook diagnostics, vectorised springs, start ramp, per-mesh calibration. 113 runs on the cluster (`~/finray_fsi/sim_data/fish-ib2d/NACA0020_Wing/`): viscosity floor 30× water, bistable deflected wing; see `RESULTS.md`.

### Validation runs of the centreline model (dt = 1e-4, 192², ~10 ms/step)

| viscosity | Re (chord) | outcome |
|---|---|---|
| 1× water (reference) | ~4e4 | NaN at t ~ 0.7 s |
| 10× | ~4e3 | NaN at t ~ 1.25 s |
| **30×** | ~1.3e3 | stable, full 10 s incl. the AoA switch |
| 100× | ~4e2 | stable, full 10 s |

| | reference | centreline 30× | centreline 100× |
|---|---|---|---|
| TE y at t = 0 | 64.4 mm | 64.4 mm | 64.4 mm |
| relaxation of the pre-bend | smooth, τ ~ 0.58 s: 30, 13, 5.5, 2 mm at 0.5, 1, 1.5, 2 s | faster, overshoots (−4 mm at 1 s) | faster (2 mm at 1 s) |
| 0° AoA, t = 1.5–5 s | TE 0.8 ± 1.2 mm, drag 0.30 N/m | −2 ± 5 mm, drag 0.8 | 7.5 ± 2 mm, drag 1.2 |
| after the switch to 15° | TE 5.2 ± 1.3 mm, Fx 0.49 N/m | 23 ± 16 mm, drag 1.4, lift 2.3 N/m | 18 ± 10 mm, drag 1.6, lift 2.1 N/m |

Qualitatively right; quantitatively the old IB2d wing moved 3–4× more than the reference (no KV damping, 30–100× lower Re).

### Unit conventions verified on the old model (still apply)
- The driver spreads `f · ds` with the **global** `ds = Lx/(2Nx)`: a physical nodal force F [N per unit depth] needs `f = F/ds`, so every spring/target constant is `K_phys/ds` (`to_ib`).
- IB2d's stock `.beam` force (python and matlab) has the end-vertex signs flipped relative to −∇(½k(cross−C)²): a net force on a kinked triple, spurious forces on a uniform-curvature arc. Not used.
