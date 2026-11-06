# NACA0020 wing: implementation decisions vs SoRoLily

Drafted 2026-10-18 from a planning discussion; reference verified at source and
stages 1-5 carried out the same day. Companion to `NOTES.md`, which holds the
model and pipeline (§ numbers below refer to it unless marked otherwise).

**Purpose.** The wing case is a test of IB2d as an FSI simulator. SoRoLily
(WaterLily_FSI + GVS: `gvs_wings/Wing_3order.jl` run by `gvs_tests.jl`, in
`~/finray_fsi/WaterLily_FSI/WaterLily_2Way_Coupling_Example/`) is validated against
experiments and is the ground truth. If IB2d, set up with the same (relative)
parameters, reproduces SoRoLily's output, IB2d is accurate enough for the
FinRay work. Where it cannot, the reason is recorded here as a **limitation
of IB2d**.

**Ground rule: the user decides every aspect of the setup.** On 2026-10-18
the user asked for the default suggestions of §5 to be adopted and the plan
to be run. Every default adopted that way is marked **[default]** in the
decisions log (§6), so it can be revisited; deviations from a default are
marked **[deviation]** with the evidence. Results: §7 and `RESULTS.md`.

---

## 1. How to read the comparison table

| status | meaning |
|---|---|
| **M**: matched | same value / same model, up to round-off or < 1 grid cell |
| **A**: approximated | IB2d can't represent it directly; the substitute and its justification are given, with the measured or expected error |
| **L**: IB2d limitation | can't be matched within the rules (§2, decision 2); accepted and reported as a finding about IB2d, with its expected effect |

**Src** = where the reference value was checked (files in
`WaterLily_2Way_Coupling_Example/`; WaterLily 1.8.0 and ParametricBodies 0.3.4
as pinned in its `Manifest.toml`, read at those versions). Every row was
checked at source on 2026-10-18; rows marked **corrected** differ from what
`NOTES.md` §3.1 said.

---

## 2. Decisions (user, 2026-10-18)

| # | decision |
|---|---|
| 1 | **The reference is the original, validated SoRoLily run** (`data_gvs/gvs_tests.h5`). SoRoLily is re-run with modified settings **only if needed**, i.e. only for a specific mismatch that can't be explained otherwise. |
| 2 | **IB2d is used through its hooks only.** The `external/IB2d` submodule stays unmodified; case-local code goes through IB2d's hooks (`user_force_model`, `please_Compute_External_Forcing`, ...). Anything that needs a solver change is a limitation (**L**). |
| 3 | **Reynolds number: as close to water as stable.** The remaining gap to water is documented as an IB2d limitation. |
| 4 | **Kelvin-Voigt damping is not reproduced** (accepted as an IB2d limitation). |
| 5 | **Solid density is not matched** (1000 vs 1064.6 kg/m³, −6 %); no mass points. |
| 6 | **Comparison outputs: all four** of TE displacement history, hydrodynamic forces, deformed shape, flow field. |
| 7 | **Reduced numerical verification**: dt_crit and an `Nx` × `mu` ladder; `stiffness_scale` calibrated per mesh; no Halton-seed campaign. |
| 8 | The SoRoLily source and `gvs_tests.h5` are on the cluster (`~/finray_fsi/WaterLily_FSI`). |
| 9 | **Proceed with the default suggestions of §5** and run the plan on the cluster (agent allowed to submit jobs for this work). |

---

## 3. Comparison table: SoRoLily → IB2d

### 3.1 Fluid and domain

| quantity | SoRoLily | IB2d | status | basis / notes | Src |
|---|---|---|---|---|---|
| fluid density | ρ = 1000 kg/m³ | `rho = 1000` | M | | `gvs_functions.jl:87` (0.5·1000·U²) |
| viscosity | ν = 1.0023e-6 m²/s (water), Re = U c/ν = 3.99e4 | **30× water** (Re 1330): the lowest viscosity that runs, at Nx 192 and 384 (`RESULTS.md` §3) | **L** (decision 3) | IB2d's advection is central-differenced (skew-symmetric, `please_Update_Fluid_Velocity.py:56`), with no upwinding: unstable at high cell Re. WaterLily's is QUICK with a median limiter (upwind-biased, dissipative), which is what lets the reference run water on 32 cells/chord | `gvs_tests.jl:40-41,57`; WaterLily `Flow.jl:4,134` |
| domain | 1.2 × 1.2 m (`sz = (192, 192)` cells, `L_c = 32` cells per chord) | `Lx = Ly = 1.2` | M | | `gvs_tests.jl:34,51` |
| grid | dx = c/32 = 6.25 mm | `Nx = 192` → dx = 6.25 mm (also 96, 384 in the ladder); Lagrangian ds = dx/2 | M | | `gvs_tests.jl:34,51` |
| boundary conditions | normal velocity = U(cos α, sin α) imposed on **all four faces, outlet included** (`exitBC = false`); tangential components zero-gradient | doubly periodic + fringe x ∈ [1.0, 1.2] m relaxing u → U(cos α, sin α) at 40 s⁻¹ + domain-mean controller 20 s⁻¹ | A | **corrected** (was "unknown"). IB2d has no inflow/outflow BCs; the fringe + mean controller is the hook-based equivalent of a uniform-flux box. Sensitivity checked in stage 5 (O3) | WaterLily `core.jl:200-219` |
| initial flow | **u = U everywhere at t = 0** (body included), then BDIM + projection | at rest, then the mean flow is driven 0 → U along a smoothstep in 0.02 s (`flow.start_ramp`) | A | **corrected** (was "unknown"). A one-step impulse (the literal equivalent) excites the clamp's ringing mode to ±4000 N/m; a 0.02 s ramp (≈ 10 ringing periods, 2 % of c/U) avoids it (O2, **[deviation]**) | WaterLily `Flow.jl:85,138` |
| advection / time stepping | QUICK + median limiter; Heun predictor-corrector; **fixed Δt = 0.01 s** (the driver overwrites WaterLily's CFL step every step; convective CFL U Δt/dx = 0.32) | central differences, FFT projection, fixed `dt` (stage 3) | **L** | **corrected** (was "adaptive CFL?"). The solver difference under evaluation | `body_simulations.jl:71`; WaterLily `Flow.jl:156-167` |
| body seen by the fluid | BDIM (ε = 1 cell); degree-1 NURBS (a polygon) through ±w/2 at every other station: 6 rigid + 10 soft per side, 31 points, 30 panels | IB points on the smooth fitted outline, spacing h = ds/2 | A | the polygon lies up to 3.7 mm (0.59 dx) inside the true outline, at the LE; its area is 1.4 % smaller. Sub-cell | `gvs_functions.jl:9-14`; `body_simulations.jl:132-140` |

### 3.2 Flow schedule

| quantity | SoRoLily | IB2d | status | basis / notes | Src |
|---|---|---|---|---|---|
| free-stream speed | U = 0.2 m/s | `flow.U = 0.2` | M | | `gvs_tests.jl:38` |
| angle of attack | 0° → 15° at t = 5 s, 0.1 s cubic smoothstep, by rotating the inflow | same schedule (`wing_case.aoa_rad`), rotating the fringe/mean target | M | | `gvs_tests.jl:44-47,72-79`; `body_simulations.jl:166-200` |
| run length | 10 s; outputs every fluid step (0.01 s), 1002 samples | `Tfinal = 10`; diagnostics every 0.01 s (hooks), fields every 0.05 s | M | | `gvs_tests.jl:64` |

### 3.3 Geometry

| quantity | SoRoLily | IB2d | status | basis / notes | Src |
|---|---|---|---|---|---|
| chord | 0.2 m (0.06 + 0.14) | `CHORD = 0.2` | M | | `Wing_3order.jl:56` |
| thickness w(x) | fitted polynomial `thickness_from_physical_x` (unclipped: −4e-5 m at the TE) | verbatim, clipped at 0 | M | | `Wing_3order.jl:65-75` |
| position | link 2 centred on the pivot: LE at 0.37, pivot (0.4, 0.6) m | same | M | | `Wing_3order.jl:98`; `gvs_tests.jl:68` |
| span | h = 0.1 m; forces × w = 0.1 | per unit depth; reference forces ÷ 0.1 to compare | M | structure and load both scale with span | `Wing_3order.jl:85`; `gvs_tests.jl:70` |

### 3.4 Structure

| quantity | SoRoLily | IB2d | status | basis / notes | Src |
|---|---|---|---|---|---|
| carriage (link 1) | `P`+`R` joints, frozen (`state_dot[7:8] .= 0`, zero initial rates) | not modelled; wing fixed | M | the h5 confirms heave = pitch = 0 throughout | `gvs_functions.jl:65` |
| rigid LE block (link 2) | 0-0.06 m, no strain DOF (its E = 1e6, η = 1e4, ρ = 625 are unused) | every point with x′ ≤ 0.06 is a target, `update_target = 0`; springs between two targets dropped | M (as k_T → ∞) | factor 3 chosen in stage 3 (O4). Clamp drift ≤ 0.2 mm; the clamp rings at ~2 ms (tether-mass mode, **IB artefact**, filtered out of the forces) | `Wing_3order.jl:28-31,146-150` |
| flexible link (link 3) | Cosserat rod, curvature only, cubic Legendre basis (κ = Σ q_k P_k(2X−1)/L₃), 4 DOF | 2D elastic solid: Halton fill + constrained Delaunay mesh, one tributary-area spring per edge | A | `stiffness_scale` fitted per mesh to the reference's own Ritz model (`verify_structure.measure_scale`, now automatic when `geom.stiffness_scale` is null). Remaining: shear near the thick root, TE curl | `Wing_3order.jl:152-160`; `SoRoLily2D.jl:584-768` |
| reference stiffness integration | K by 5-point Gauss quadrature | calibration target = exact integral | M | Gauss vs exact: tip deflection −0.14 %, at 0.9 L₃ −0.02 % | `Wing_3order.jl:41`; `SoRoLily2D.jl:506,906-958` |
| Young's modulus | E = 1e4 Pa | `E_flex = 1e4` × measured `stiffness_scale` (4.083 at ×4, Nx 192) | A | | `Wing_3order.jl:28` |
| Poisson's ratio | ν = 0.5 (only enters shear, which is frozen) | effective 0.29-0.34 | **L** | central-spring network limit | `Wing_3order.jl:29` |
| material damping | Kelvin-Voigt η = 2e3 Pa·s; bending damping 3η I → τ = 3η/E = 0.6 s | none | **L** (decision 4) | expected: faster relaxation, overshoot, larger TE fluctuations | `Wing_3order.jl:30`; `SoRoLily2D.jl:868-876` |
| solid density | ρ_s = 1064.6 kg/m³ (fluid excluded from the body by BDIM) | interior fluid at 1000, moving with the wing | **L** (decision 5) | −6 % inertia | `Wing_3order.jl:31` |
| actuation | tendons, tension × 0 | none | M | | `gvs_functions.jl:63` |
| initial state | q = [1, 0, 0, 0]: κ = 1/L₃, TE 64.36 mm up, at rest, released at t = 0 | bent `.vertex`, straight rest lengths: TE 64.36 mm | M | | `gvs_tests.jl:60`; h5 |

### 3.5 Coupling and numerics

| quantity | SoRoLily | IB2d | status | basis / notes | Src |
|---|---|---|---|---|---|
| coupling scheme | **partitioned, explicit (staggered)**: per fluid step, panel pressures → 6th-order Butterworth low-pass at **3 Hz**; panel forces (pressure + viscous) → low-pass at **1.5 Hz**; GVS ODE integrated over the step (Rodas5) with that load held fixed; body updated | IB: massless structure, unfiltered forces spread every step | **L** | **corrected** (was "to check"). The reference's structure only feels loads below ~1.5 Hz, delayed by the causal filters (forces start ~0.3 s late, visible in `reference/reference_overview.png`). IB2d can't filter its own coupling | `body_simulations.jl:34-46,60-100,144-150` |
| time step | fixed 0.01 s | fixed `dt` from stage 3 | — | | `gvs_tests.jl:29` |
| mesh density | n/a | multiplier ×4 at every `Nx` (h/dx fixed) | **[default]** O4 | calibration within 2.1 % of ×32 | — |
| spring force cost | — | IB2d's own spring force is a Python loop (17 µs/spring/step on the cluster); computed instead in the user-force hook, vectorised, same formula (`spring_force_check.py`: agrees to 2e-16) | **[default]** O6 | step time 157 → ~60-80 ms at Nx 192 | — |

### 3.6 Outputs (decision 6)

| output | SoRoLily | IB2d | status | notes | Src |
|---|---|---|---|---|---|
| TE displacement y(t) | TE = last soft backbone station, from q(t) via the reference's own kinematics (`reference/sorolily.py`: midpoint rule, 20 segments) | TE vertex, from the straight wing | M | windows: release 0-2 s, 0° 2-5 s, 15° 5.5-10 s | `SoRoLily2D.jl:1556-1630` |
| forces | `force_history` rows: [Σ panel torques about their own reference points, **Σ Fx**, 6 generalized forces]; **total Fy is not stored** | tether forces on the clamp (whole wing, pressure + viscous), every step, then the reference's own filter (Butterworth 6, 1.5 Hz) | A | **corrected**: Fy rebuilt from the stored panel pressures (`reference/extract_reference.py`); the same rebuild of Fx matches the stored Fx to 0.3 mN RMS, so it is exact and the viscous part is negligible | `gvs_functions.jl:70-100`; `gvs_tests.jl:93-97` |
| deformed shape | centreline at the 21 soft stations from q(t) | same material points (rest centreline mapped through the deformed triangulation) | M | | — |
| flow field | **not stored**: only the filtered pressure at the 30 panel midpoints (gauge: the corner cell, `p .-= p[2,2]`) | surface pressure at the same panel midpoints, same gauge point, every field frame; vorticity snapshots (IB2d only) | A | **[default]** O8: compared through surface pressure; whole-field comparison would need a SoRoLily re-run (decision 1: not triggered) | `body_simulations.jl:74,77-78`; `gvs_tests.jl:96` |

---

## 4. Work carried out (2026-10-18)

| stage | what | where |
|---|---|---|
| **0. Verify the reference at source** | every row of §3 (corrections marked) | this file |
| **1. Reconcile the setup** | `stiffness_scale: null` → measured for the exact mesh at geometry time (4.083 at ×4, Nx 192; the old 4 / 3.966 were for other meshes); multiplier ×4; target factor and dt from stage 3 | `Wing_Geom.py`, `verify_structure.measure_scale` |
| **2. Tooling** | reference extraction (`reference/`), per-step diagnostics through the hooks (`wing_diagnostics.py`), start-up ramp, vectorised springs, `postprocess.py` (all four outputs vs reference), `sweeps.py` + `slurm/` (as FinRay), `report.py` | case folder |
| **3. dt_crit** | 1 s runs on a dt ladder per `Nx`, target factor 3 and 15, mu 30× (and 1× at Nx 192) | sweep `stage3_dt` |
| **4. `Nx` × `mu` ladder** | full 10 s schedule; then, because everything below 30× diverged: diagnosis (halved dt, no pre-bend), refinement of the floor (15×, 20×), structural mesh ×1 | sweeps `stage4_ladder`, `stage4_diag`, `stage4_refine`, `stage4_mesh` |
| **5. Comparison runs** | Nx 192 and 384 at 30× (from the ladder) + sensitivity: fringe ×0.5/×2 (O3), no ramp (O2), factor 15 (O4), mesh ×8, 100× water, no pre-bend, pre-bend −1 / +0.5 | sweeps `stage5_sensitivity`, `stage5_bend_sign` |
| 6. SoRoLily re-run | not triggered | — |

Results are in `~/finray_fsi/sim_data/fish-ib2d/NACA0020_Wing/<sweep>/` (each with
`report/summary.md`) and summarised in `RESULTS.md` (this folder).

---

## 5. Open decisions and the defaults adopted

| # | decision | default (adopted 2026-10-18) |
|---|---|---|
| O1 | metric definitions and tolerances | all measures computed (`postprocess.py`); no pass/fail tolerance for the IB2d-vs-SoRoLily comparison (still the user's call); numerical convergence judged at 2 % (FinRay's default) |
| O2 | initial flow state | match the reference's start (u = U at t = 0); implemented as a 0.02 s ramp, not a one-step impulse (**[deviation]**, §6) |
| O3 | fringe parameters | keep 40 / 20 s⁻¹; one sensitivity pair (×0.5, ×2) at the chosen setting |
| O4 | baseline mesh and clamp | ×4; target factor from stage 3 |
| O5 | `mu`, `Nx` levels; production dt | mu ∈ {1, 3, 10, 30} × water, Nx ∈ {96, 192, 384} (**[deviation]**: 96 replaces 768, which would exceed the 2-day limit); dt = largest ladder value ≤ 0.7 dt_crit |
| O6 | vectorised spring force | yes (the trigger, stage-4 runs > 1 day, was met on the cluster's nodes) |
| O7 | reading `gvs_tests.h5` | `h5py` in the `ib2d` env; the extract is committed as `reference/sorolily_reference.npz` |
| O8 | flow-field comparison | surface pressure at the reference's panels (all the h5 holds); no SoRoLily re-run |
| O9 | mesh noise | accepted (no seed runs) |

---

## 6. Decisions log

| date | # | decision | reason |
|---|---|---|---|
| 2026-10-18 | 1-8 | see §2 | planning discussion |
| 2026-10-18 | 9 | proceed with the §5 defaults; run on the cluster | user |
| 2026-10-18 | stage 0 | §3 rows verified; corrections: BCs, initial flow, fixed 0.01 s step, QUICK advection, filtered staggered coupling, force_history layout (no Fy), pressure gauge | source read at the pinned versions |
| 2026-10-18 | O1 **[default]** | no pass/fail tolerance for IB2d vs SoRoLily; numerical convergence at 2 % | tolerance for the physics comparison stays the user's decision |
| 2026-10-18 | O2 **[deviation]** | 0.02 s smoothstep ramp of the mean flow instead of a one-step impulse | the impulse excites the clamp's ringing mode (period ~1.8 ms) to ±4000 N/m, undamped over 0.1 s (`_dev` runs); with the ramp the ringing is 29 N/m RMS, as without any forced start (21 N/m), and filtered out of the force comparison |
| 2026-10-18 | O4 **[default]** | multiplier ×4; `stiffness_scale` measured per mesh | ×4 within 2.1 % of ×32 |
| 2026-10-18 | O5 **[deviation]** | Nx ladder {96, 192, 384} | 768 would exceed the 2-day limit; 96 gives the third grid for a convergence estimate |
| 2026-10-18 | O6 **[default]** | springs computed vectorised in the user-force hook | 17 µs/spring/step on the cluster: Nx 384 runs projected beyond 1 day |
| 2026-10-18 | O7, O8, O9 **[default]** | see §5 | |
| 2026-10-18 | O4 (stage 3) | **target factor 3**: dt_crit 1e-3 / 3.125e-4 / 1.5625e-4 s at Nx 96 / 192 / 384; factor 15 halves it; production dt 6.25e-4 / 2e-4 / 1e-4 s (≤ 0.7 dt_crit) | larger dt at every grid; the clamp's static effect is 0.04 %, its dynamic effect +15 % on the 0° offset (`stage5_sensitivity`) |
| 2026-10-18 | O6 | vectorised springs adopted (`spring_force_check.py`: 2e-16) | step 157 → ~24 ms at Nx 192 on the batch nodes |
| 2026-10-18 | O5 / decision 3 | **viscosity: 30× water** (Re 1330), the lowest that completes the schedule (Nx 192, 384); 20× (Nx 192), 15× and 20× (Nx 384), and every level ≤ 10× diverge | grid-scale checkerboard of IB2d's collocated central-difference solver; not dt, not the release, not the structural mesh (`stage4_diag`, `stage4_mesh`) |
| 2026-10-18 | additions | stages 4-5 extended beyond the plan (`stage4_diag`, `stage4_refine`, `stage4_mesh`, extra `stage5_sensitivity` variants, `stage5_bend_sign`) | the ladder and the first comparison showed results the plan did not anticipate (viscosity floor, a bistable deflected wing); each extension is one targeted test |
| 2026-10-18 | O8 (pressure sampling) | surface pressure sampled supp/2 · dx outside the IB outline | on the outline IB2d's pressure mixes the interior (±800 Pa) |
| 2026-10-18 | decision 1 | **trigger for a SoRoLily re-run is met but the re-run was not done** (recommended in `RESULTS.md` §5: ν × 30, and ν × 30 with η = 0) | changing the ground-truth run and installing its toolchain is the user's call |
| 2026-10-20 | decision 1 (user) | **SoRoLily re-run at IB2d's working Re** (ν × 30), with and without its damping; η = 0 diverges, so η was swept (1000, 500, 200, 100) and η = 200 (lowest stable) stands in for "undamped" | user request; results in `RESULTS.md` §4b and the report: SoRoLily stays straight at matched Re, so IB2d's deflection is an IB2d artefact |
| 2026-10-18 | tooling | post-processing also runs on crashed runs; IB2d stops with an IndexError on NaN positions instead of writing NaN frames, so stability is judged from the per-step tether force (`analysis.stable`) | |

---

## 7. Results

See `RESULTS.md`.
