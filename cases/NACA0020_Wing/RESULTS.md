# NACA0020 wing: IB2d vs SoRoLily — results

Runs of 2026-10-18 on the ICF cluster, following `IMPLEMENTATION_DECISIONS.md`
(§4 stages, §5 defaults). Data, per-run figures and per-sweep reports:
`~/finray_fsi/sim_data/fish-ib2d/NACA0020_Wing/` (its `README.md` explains the layout;
cross-sweep figures in `SUMMARY/`). Every number below comes from the runs'
`run_info.json` (`analysis`, computed by `postprocess.py`), or from
`SUMMARY/runs.csv`.

## 1. Verdict

**At the resolutions it can afford, IB2d does not reproduce the SoRoLily wing,
and it cannot run this case anywhere near the reference's Reynolds number.**

1. **Viscosity floor.** Stock IB2d (hooks only) runs the 10 s schedule only at
   **30× water's viscosity or more** (Re 1330 instead of 4.0e4), at Nx = 192
   (32 cells/chord, the reference's grid) and Nx = 384. Every run below that
   diverges within 0.2-2 s, at every grid tested (96-384), every dt, both
   structural meshes, with or without the pre-bend's release (§3).
2. **Marginal even there.** At 30× water and Nx 192, two of four initial
   conditions diverge (straight start at 0.51 s, half pre-bend at 1.14 s), and
   the coarser structural mesh (×1) diverges at 1.4 s. Only the full ±pre-bend
   runs complete.
3. **Where it runs, the wing behaves differently** (§4): the release is ~3×
   faster and overshoots (no Kelvin-Voigt damping), and afterwards the wing
   locks into a **deflected state, +22 mm at the trailing edge at 0°** (the
   reference: 0.3 mm), whose sign follows the pre-bend's (−24 mm with a
   downward pre-bend) and which survives the switch to 15° — a bistable
   (statically divergent) wing. It persists under grid refinement (Nx 384:
   +29 mm at 0°) and is insensitive to the fringe, the start-up, the
   structural mesh, the clamp stiffness (within 15 %) and Re between 400 and
   1330.
4. **Most likely causes** (not yet separated, see §5): (a) the viscosity
   floor — at Re ≤ 1330 the flow over a thick section and its flexible aft
   part is qualitatively different (the fluid pushes the deflected flap
   further, +0.5 N/m on the flap against −3.6 N/m on the block); (b) the
   missing material damping and the reference's 1.5 Hz load filtering, which
   make the reference's wing overdamped. Both are decision-3/4 limitations of
   IB2d (no upwinding; no implicit damping).

## 2. What was run

| sweep | runs | purpose | outcome |
|---|---|---|---|
| `stage3_dt` | 76 × 1 s | dt_crit per Nx, target factor 3 vs 15 | dt_crit (factor 3) = 1e-3, 3.125e-4, 1.5625e-4 s at Nx 96, 192, 384; factor 15 halves it → **factor 3** (static effect of the clamp: 0.04 %) |
| `stage4_ladder` | 12 × 10 s | Nx {96, 192, 384} × mu {1, 3, 10, 30} × water | stable: Nx 192 and Nx 384 at 30×; all others diverge (§3) |
| `stage4_diag` | 5 × 2 s | halved dt; no pre-bend | divergence time unchanged by dt (0.60 vs 0.65 s; 0.2135 vs 0.2134 s): **not a time-step limit**; earlier without pre-bend |
| `stage4_refine` | 3 × 10 s | 20× at Nx 192, 15× / 20× at Nx 384 | all diverge (0.94, 0.26, 0.50 s) |
| `stage4_mesh` | 8 × 10 s | structural mesh ×1 (Lagrangian spacing h = IB2d's ds) | all diverge, even 30× (1.4 s): the dense ×4 mesh is not the cause |
| `stage5_sensitivity` | 7 × 10 s | at Nx 192, 30×: fringe ×0.5/×2, no pre-bend, 100×, no ramp, factor 15, mesh ×8 | §4 table |
| `stage5_bend_sign` | 2 × 10 s | pre-bend −1 and +0.5 | −1: offset flips (−24 mm); +0.5: diverges at 1.14 s |

Production settings (`case_params.json` + `sweeps.DT_PRODUCTION`): mesh ×4,
`stiffness_scale` measured per mesh (4.083 at Nx 192, 4.055 at Nx 384), target
factor 3, dt = 6.25e-4 / 2e-4 / 1e-4 s at Nx 96 / 192 / 384 (largest ladder
value ≤ 0.7 dt_crit), vectorised springs in the hook, start ramp 0.02 s.
Cost on the cluster's nodes: 24 ms/step at Nx 192 (20 min per 10 s run),
~150 ms/step at Nx 384 (~4 h); peak memory < 0.8 GB per run.

## 3. The viscosity floor

Figure: `~/finray_fsi/sim_data/fish-ib2d/NACA0020_Wing/SUMMARY/stability_map.png` (every run in the Nx × viscosity plane).

- Below 30× water the runs blow up through a **grid-scale (odd-even)
  checkerboard** in pressure and vorticity that starts in the IB layer around
  the wing and radiates into the domain. IB2d's fluid solver is collocated
  with central differences (`please_Update_Fluid_Velocity.py`): its projection
  cannot see the Nyquist modes, nothing in the advection damps them, and only
  viscosity does. WaterLily (the reference) uses QUICK with a limiter
  (upwind-biased, dissipative), which is what lets it run water on the same
  32 cells/chord.
- It is **not a time-step limit** (halving dt changes nothing), **not the
  release** (a straight start diverges earlier), **not the structural mesh**
  (×1 is worse, ×8 the same).
- In cell Reynolds number U·dx/ν: stable at 42 (Nx 192, 30×) and 13 (Nx 192,
  100×); diverged at 62, 83, 125 and above — but also at 42 on Nx 384 (15×),
  and finer grids blow up *sooner* (the grid-scale modes grow at a rate
  ~U/dx). So it is not a single cell-Re threshold; empirically the stable
  region is **≥ 30× water at Nx 192-384** (Nx 96 fails even there, at 2 s).
  At Nx 384 / 30× the grid-scale KE still grows slowly (+0.04 s⁻¹ over the
  run; −0.14 s⁻¹ at Nx 192): even that setting is close to the edge.
- What water would need: if the threshold were a cell Re of ~40, dx ≈ 0.2 mm,
  ~1000 cells per chord (Nx ≈ 6000) — far beyond IB2d's cost (a 10 s run at
  Nx 384 already takes 4 h); and since 15× fails at Nx 384, even that may not
  suffice.

## 4. The comparison at the stable setting (Nx 192, 30× water, dt 2e-4)

Windows: release 0-2 s, 0° 2-5 s, 15° 5.5-10 s. Forces per unit depth,
1.5 Hz low-pass (the reference's own filter). SoRoLily values from
`reference/sorolily_reference.npz` (`SUMMARY/reference_metrics.json`).

| measure | SoRoLily (Re 4.0e4) | **IB2d Nx 192** | **IB2d Nx 384** | fringe ×0.5 | fringe ×2 | no ramp | factor 15 | mesh ×8 | 100× (Re 400) | pre-bend −1 |
|---|---|---|---|---|---|---|---|---|---|---|
| status | — | stable | stable | stable | stable | stable | stable | stable | stable | stable |
| TE t50 (release) [s] | 0.46 | 0.16 | 0.14 | 0.16 | 0.16 | 0.20 | 0.16 | 0.16 | 0.16 | n/a |
| TE minimum, release [mm] | +2.1 | −5.2 | +7.1 | −5.2 | −5.3 | +1.9 | −3.3 | −6.5 | +9.4 | n/a |
| TE mean ± std, 0° [mm] | 0.3 ± 0.5 | 21.8 ± 6.0 | 29.4 ± 3.4 | 21.8 ± 7.1 | 22.2 ± 4.9 | 21.4 ± 12.3 | 25.2 ± 5.4 | 21.8 ± 8.7 | 26.6 ± 5.5 | −24.2 ± 8.2 |
| TE mean ± std, 15° [mm] | 4.9 ± 1.7 | 38.8 ± 6.2 | 37.1 ± 5.9 | 35.2 ± 9.9 | 36.9 ± 11.7 | 31.1 ± 13.0 | 36.5 ± 13.0 | 37.4 ± 8.8 | 30.0 ± 7.1 | −21.2 ± 9.5 |
| Fx mean, 0° [N/m] | 0.30 | 1.36 | 1.07 | 1.39 | 1.35 | 1.60 | 1.20 | 1.42 | 1.72 | 1.52 |
| Fy mean, 0° [N/m] | 0.54 | −3.03 | −1.63 | −3.05 | −3.02 | −2.42 | −3.31 | −2.75 | −1.67 | +2.16 |
| Fx mean, 15° [N/m] | 0.45 | 0.72 | 0.32 | 0.96 | 0.95 | 1.23 | 0.72 | 0.78 | 1.18 | 1.68 |
| Fy mean, 15° [N/m] | 2.37 | 3.22 | 2.80 | 2.39 | 2.61 | 1.68 | 3.61 | 3.12 | 2.01 | 4.30 |
| centreline RMS distance, t = 0.5 s [mm] | — | 12.6 | 9.0 | 12.6 | 12.6 | 9.7 | 12.4 | 13.1 | 8.2 | 16.4 |
| centreline RMS distance, mean shape 0° / 15° [mm] | — | 11.6 / 16.6 | 14.5 / 16.2 | — / 15.5 | — / 16.1 | — / 13.5 | — / 16.3 | — / 15.9 | — / 13.1 | 13.2 / 13.3 |
| surface p, RMS diff. of the window mean 0° / 15° [Pa] (reference RMS 9.6 / 13.5) | — | 10.1 / 5.5 | 7.6 / 6.9 | — / 5.3 | — / 5.0 | — / 5.5 | — / 5.5 | — / 4.9 | — / 5.5 | 5.7 / 5.6 |
| ms / step (cluster) | — | 24.5 | 145 | 22.0 | 21.9 | 24.1 | 22.1 | 30.9 | 24.1 | 20.2 |

(—: see each sweep's `report/summary.md` for every measure of every run; t50 and the release minimum are not defined for a downward pre-bend.)

- **Release.** The undamped IB2d wing straightens ~3× faster (t50 0.14-0.16 s
  vs 0.46 s), overshoots to −5 mm at 0.55 s (Nx 192; Nx 384 dips to +7 mm) and
  rebounds to +30 mm by 1 s; the
  reference relaxes monotonically (Kelvin-Voigt τ = 0.6 s, and loads filtered
  at 1.5 Hz, so its structure feels the flow ~0.3 s late).
- **Deflected state.** After the release, the IB2d wing settles at ±22 mm on
  the side the release leaves it, bent smoothly along the whole flexible link
  (1, 6, 14, 22 mm at s/L₃ = 0.14, 0.43, 0.71, 1) — not a trailing-edge
  curl, no inverted triangles. The flexible part carries only +0.5 N/m of net
  fluid load, the block −3.6 N/m: the deflection is held by the fluid's
  moment on the flap, a static divergence absent from the reference.
- **Forces.** Drag at 0° is 3.6-4.5× the reference's (viscous drag at Re 1330
  and the deflected shape); lift at 15° is the closest measure (+18 % at
  Nx 384, +36 % at Nx 192), but at 0° the sign is wrong (the deflected wing is
  negatively cambered).
- **Grid dependence (Nx 192 → 384, 30×).** 15° measures move by 4-13 % (TE
  mean, Fy), the 0° offset by +35 % and the 0° lift by a factor 2: the
  comparison is not grid-converged, but the qualitative picture (deflected,
  bistable wing; fast undamped release) is the same on both grids.
- **Flow.** The reference stored no flow fields; the surface pressure at its
  30 panels is the only flow quantity to compare (§3.6 of the decisions). IB2d's
  vorticity (`flow_snapshots.png` in each run) shows a separated, shedding wake
  behind the deflected flap.

## 4b. Matched-Re check (2026-10-20)

SoRoLily was re-run at IB2d's working Re (ν × 30, Re 1330) with its own code
(installed on the cluster; the water re-run reproduces the stored reference to
2 µm), and with its wing damping reduced from 2000 to 1000, 500, 200, 100 and
0 Pa·s. Runs and README: `~/finray_fsi/sim_data/fish-ib2d/NACA0020_Wing/sorolily_matched_re/`.

- At Re 1330, SoRoLily's wing stays straight: TE 0.3 mm at 0°, 3.6 mm at 15°
  (water: 0.3 / 4.9). With η = 200 the release is faster than IB2d's
  (t50 0.07 s) and overshoots as far (−4.8 mm), yet it returns to neutral
  (0.3 / 4.1 mm). **IB2d's 22-29 mm deflection is therefore an IB2d artefact**,
  not low-Re physics or the missing damping.
- Matched-Re forces, SoRoLily (η 2000 / 200) vs IB2d Nx 192 / 384: drag 0°
  0.40 / 0.40 vs 1.36 / 1.07 N/m; lift 15° 1.86 / 1.83 vs 3.22 / 2.80 N/m.
- SoRoLily diverges at η = 100 (0.2 s) and η = 0 (0.11 s); the η = 0 wing
  without fluid is stable, so its explicit staggered coupling needs the
  damping (light wing, ρ_s/ρ ≈ 1.06). Part of η = 2000 may be numerical.
- IB2d comparisons: `comparisons.sorolily_nu30[_eta200]` in the Nx192/384
  30× runs' run_info.json.

## 5. Recommended next steps (for the user to decide)

0. *(Done 2026-10-20, §4b: the deflected state is an IB2d artefact. Next:
   a rigid wing in both codes at Re 1330, to separate IB2d's flow solution
   from its fluid-structure model.)*
1. **Separate the two causes with a SoRoLily re-run** (decision 1's trigger is
   met: the deflected state can't be attributed from IB2d alone). Same
   setup with (a) ν × 30, and (b) ν × 30 and η = 0. If (a) also diverges into
   a deflected state, the cause is the Re floor and IB2d's FSI coupling is not
   at fault for it; if (b) does and (a) doesn't, it is the missing damping.
   Needs Julia + the pinned packages (`Manifest.toml`, Julia 1.13); the run is
   1000 fluid steps on 192², cheap. (With ν × 30 the reference's fixed 0.01 s
   step stays within WaterLily's CFL limit.)
2. **Decide the tolerance (O1)**: the gaps above (TE offsets of 20-40 mm on a
   0.2 m chord, forces off by 35-450 %) exceed any reasonable tolerance, so
   for this case the answer does not depend on it.
3. For FinRay, the same floor applies: plan for ≥ 30× water (or a solver
   change: upwinding/filtering in IB2d's advection, which decision 2 rules
   out for the evaluation of stock IB2d).

## 6. Caveats

- Defaults adopted without the user's case-by-case review are marked
  **[default]** in `IMPLEMENTATION_DECISIONS.md` §6; deviations **[deviation]**
  (the start ramp; Nx 96 instead of 768). Stages 4-5 were extended beyond the
  plan (`stage4_diag`, `stage4_refine`, `stage4_mesh`, the extra stage-5
  variants, `stage5_bend_sign`) once the ladder showed results the plan did not
  anticipate.
- Grid convergence of the comparison measures is not established: only Nx 192
  and 384 run at 30× (two grids, no order or GCI), and they differ by up to
  35 % (0° offset). No Halton-seed spread (O9).
- dt_crit came from 1 s runs; all production runs that diverged did so for
  the viscosity reason above, not dt (stage4_diag).
- The IB2d forces are the clamp tether forces: they include the clamp's
  ringing (~2 ms period, 6-30 N/m RMS), filtered out at 1.5 Hz like the
  reference's forces.
- Surface pressure is sampled 2 dx (the delta support) outside the IB
  outline; on the outline itself IB2d's pressure mixes in the interior
  (±800 Pa).
