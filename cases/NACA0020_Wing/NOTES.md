# NACA0020 wing – WaterLily_FSI (GVS) → IB2d

Replicates `WaterLily_2Way_Coupling_Example/gvs_wings/Wing_3order.jl` as run by `gvs_tests.jl`.

## Run

```sh
conda activate ib2d          # numpy, numba, pyfftw, vtk, matplotlib
cd pyIB2d/Examples/NACA0020_Wing
python Wing_Geom.py --plot   # writes wing.vertex/.spring/.user_force/.target/.geo_connect + geometry_preview.png
python verify_structure.py   # fluid-free calibration check (see below)
python main2d.py             # the FSI run  (~10 ms/step at 192x192; 10 s at dt=1e-4 is ~17 min)
python analyze_run.py        # tip motion + lift/drag history from viz_IB2d
```
All tunables are in `case_params.json`; grid/fluid/time step are read from `input2d` (single source of truth for `Wing_Geom.py`).

## What the reference model actually is

| GVS (`Wing_3order.jl` + `gvs_functions.jl`) | Consequence |
|---|---|
| link 1 = massless carriage, `P`+`R` joints, but `dynamics_gvs` sets `state_dot[7:8] = 0` and `init_state` has zero velocity | heave and pitch are **frozen**: the wing is a fixed body |
| link 2, 0.06 m, `Phi_dof = 0` → rigid | fixed leading-edge block, 0 – 0.06 m from the LE (30 % chord) |
| link 3, 0.14 m, curvature only, cubic Legendre basis (`Phi_dof=(1,0,0)`, `Phi_odr=(3,0,0)`) | **inextensible, unshearable Euler–Bernoulli cantilever**, clamped at 0.06 m; stretch/shear stiffness (`E A`, `G A`) never enters |
| `E = 1e4 Pa`, `ρ = 1064.6`, `η = 2e3 Pa s`, `ν = 0.5` (`G` unused) | EI(x) = E·h·w(x)³/12; KV bending damping 3·η·I (not reproduced, see limitations) |
| `w(x)` fitted thickness, zero at the tip, ≈ 0.04 m near x = 0.06 m | NACA0020, chord 0.2 m |
| `h = 0.1 m` span, forces × `w = 0.1` | structure and fluid load both ∝ span ⇒ response independent of span ⇒ **IB2d uses span = 1 (per unit depth)** |
| tendons `n_sact = 2` | `dynamics_gvs` multiplies the tension by 0 → **no actuation** (not modelled) |
| `init_state = [0,0,1,0,0,0, 0…]` | starts **pre-bent**: constant-curvature mode = 1 ⇒ κ = 1/L₃ = 7.14 m⁻¹, tip rotated 1 rad (57°), released at t = 0 (`init_q` in the json) |
| inflow `U = 0.2 m/s`, α: 0° → 15° at t = 5 s (0.1 s smoothstep), ν = 1.0023e-6, ρ = 1000; wing pivot at (0.4, 0.6) m in a 1.2 × 1.2 m box; 32 cells / chord | reproduced (fringe forcing, 192 × 192 grid, same pivot) |

## Mapping to IB2d fibres

| Quantity | IB2d representation |
|---|---|
| airfoil outline | **skin**: 2 rails of vertices (upper/lower), arc-length spacing ds = Lx/2Nx, LE/TE closed by the centerline end vertices |
| rigid LE block + clamp | **target points** (`K_target`) on every vertex with x′ ≤ 0.06 m |
| EB bending, EI(x) | coarse **centerline** chain (spacing `centerline_spacing`, 10 mm) with a turning-angle energy ½(EI/s)θ² per vertex (`.user_force`); EI is the harmonic mean of E·w³/12 over the vertex's cell; the clamp vertex gets 2× (its rotation represents half a cell) |
| inextensibility | axial **springs** on the centerline, k = E·w/s |
| skin follows the section (plane sections stay plane, no shear) | each skin vertex is tied to both ends of its centerline segment by 2 stiff **springs** (rigid triangle) |
| outline draw / even spacing | very weak outline springs (`K_skin`, changes EI by < 1 %) |
| initial pre-bend | vertices start bent; rest lengths / rest angles are those of the straight wing, so the bend is the initial elastic state |
| free stream, AoA(t) | `please_Compute_External_Forcing.py`: fringe (sponge) forcing at x ∈ [1.0, 1.2] relaxing u → U(cos α, sin α) at 40 s⁻¹ (IB2d is doubly periodic; there is no inflow BC) |

### IB2d unit conventions that matter (verified in `IBM_Driver.py`)
* The driver spreads `f · ds` with the **global** `ds = Lx/(2Nx)`. A physical nodal force F [N per unit depth] therefore needs `f = F/ds`; every spring/target/bending constant is `K_phys / ds` (`Wing_Geom.to_ib`, `bend_k`). Note `FinRay_Geom.py` uses `k = E·A/ds` without the extra `1/ds`, i.e. its stiffnesses are effectively ds× larger than the material values.
* **IB2d's stock `.beam` force is not conservative** (python and matlab): the left/right-vertex signs are flipped relative to `-∇(½k(cross-C)²)`, so a triple carries a net force (measured: kinked triple → net (−0.025, −1.0), and a uniform-curvature arc, i.e. a pure moment, feels a spurious force on every interior vertex). `k_beam` can therefore not be calibrated against EI, and bending is done through `user_force_model` instead (`give_Me_General_User_Defined_Force_Densities.py`; conservative, zero net force/moment, checked by finite differences). This also affects the small "regularisation" beams in `FinRay`; and the `k_Beam = E·I/ds` formula noted in `FinRay/NOTES.md` should be `EI/(ds·s⁵)` for IB2d's cross-product energy if the stock beam is ever used.

## Calibration (`verify_structure.py`, fluid-free, IB2d's own force routines)

Uniform load 0.5 N/m on the flexible link, tip deflection [mm]:

| model | tip |
|---|---|
| Euler–Bernoulli, same EI(x) | 0.770 |
| Ritz with the GVS cubic-curvature basis (what `Wing_3order.jl` solves) | 0.757 |
| IB2d chain, load on centerline (10 mm) | 0.772 (+0.2 % vs EB) |
| IB2d chain, load on the upper skin (through the ties) | 0.74 (−4 %) |

## Results of the validation runs (dt = 1e-4, 192 x 192, ~10 ms/step, 10 s = ~20 min)

| viscosity | Re (chord) | outcome |
|---|---|---|
| 1 x water (reference) | ~4e4 | NaN at t ~ 0.7 s (fluid solver; cell Re ~ 1250, central differencing) |
| 10 x | ~4e3 | NaN at t ~ 1.25 s |
| **30 x (default)** | ~1.3e3 | stable, full 10 s incl. the 0 -> 15 deg switch at t = 5 s |
| 100 x | ~4e2 | stable, full 10 s (smoothest; use as the fallback) |

Comparison with the reference run (`WaterLily_2Way_Coupling_Example/data_gvs/gvs_tests.h5`, tip deflection rebuilt from its stored curvature modes, x-force / span):

| | reference | IB2d 30 x | IB2d 100 x |
|---|---|---|---|
| tip y at t = 0 | 64.4 mm | 64.4 mm (geometry / initial state match) | 64.4 mm |
| relaxation of the pre-bend | smooth, tau ~ 0.58 s (= 3 eta / E, KV-dominated): 30, 13, 5.5, 2 mm at 0.5, 1, 1.5, 2 s | faster, overshoots (-4 mm at 1 s) | faster (2 mm at 1 s) |
| 0 deg AoA, t = 1.5 - 5 s | tip 0.8 +/- 1.2 mm, drag 0.30 N/m | -2 +/- 5 mm (wake unsteady), drag 0.8 | 7.5 +/- 2 mm, drag 1.2 |
| after the switch to 15 deg | tip 5.2 +/- 1.3 mm, Fx 0.49 N/m | 23 +/- 16 mm, drag 1.4, lift 2.3 N/m | 18 +/- 10 mm, drag 1.6, lift 2.1 N/m |

Qualitatively right (same start, deflection and positive lift after the switch, same order of magnitude); quantitatively the IB2d wing moves 3-4x more than the reference. Two causes, both understood: no Kelvin-Voigt damping (below) and ~30-100x lower Re (higher drag, unsteady wake at 0 deg).

## Deliberate deviations / known limitations
* **Kelvin-Voigt damping is NOT reproduced.** The reference is overdamped (3 eta / E = 0.6 s), so its wing is a quasi-static low-pass filter of the flow loads. In stock IB2d the structure is massless and forces are explicit, so a dashpot has loop gain c dt g / m ~ 1e2 - 1e3 (nodal inertia ~ rho dx^2). Tried and rejected: a per-vertex dashpot with a stability limiter (realises ~0.1 - 20 % of eta) and a modal dashpot in the reference's 4 cubic modes (NaN at once: the modal moment ends up as forces on the few end vertices). Reproducing it needs implicit fluid-structure coupling or a partitioned scheme (reduced GVS ODE driving target-point positions from the fluid force, as WaterLily does) - which departs from the fibre-model framework, so it is left as an open choice.
* **Reynolds number.** See table: physical-water Re is unreachable with IB2d's central-difference solver at this resolution (WaterLily uses upwind QUICK). `mu` in `input2d` is the knob. Dynamic-pressure-based quantities (loads on the structure, Cauchy number) are unchanged; viscous drag and the wake are not.
* **EI floor.** True EI ~ w^3 -> 0 at the trailing edge, so uncorrected Euler-Bernoulli has a log-divergent tip curvature under pressure (EB tip 4.9 mm vs 0.77 mm floored); the reference's cubic curvature basis regularises this implicitly. `EI_floor_thickness` (8 mm) clamps the thickness used for EI/EA near the tip; the GVS-cubic result changes < 5 % with or without it.
* **Density.** The airfoil interior is fluid (rho = 1000) that moves with the wing vs rho_s = 1064.6 in the reference (+6 %): effectively neutrally buoyant. No mass points.
* **Free stream.** Fringe forcing plus a domain-mean-velocity controller (`mean_rate`), not an inflow BC; the domain is periodic in y too. Without the mean controller the turned flow reached the wing ~2 s after the AoA switch; with it, within ~0.2 s. The wing's own upwash raises the measured upstream angle to ~18 deg after the switch.
* **Time step.** All validation used dt = 1e-4; dt = 2e-4 was only tried at unstable Re, so it is untested where the run is otherwise stable.
* **Tendons** (`n_sact = 2`) are zeroed in the reference's coupled run and are not modelled.
* Lift/drag: `analyze_run.py` uses the clamp tether reaction per unit depth; multiply by 0.1 m for a WaterLily-style force. NB the reference's `force_history[2]` (plotted as "Lift") is actually the summed x-force.
