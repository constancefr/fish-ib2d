"""
Fluid-free stiffness calibration of the FinRay wall (stage B of
EXPERIMENT_PLAN.md): measures `ribbon_fill_stiffness_scale`, the single
multiplier on every spring, so that a meshed wall bends like the real
material. Modelled on cases/NACA0020_Wing/verify_structure.py.

Why it is needed: tributary-area springs (mesh_utils.tributary_springs) put
each triangle's area into the springs exactly once, but a spring network at
scale 1 is only ~1/3 as stiff as a solid of modulus E, because each spring
resists stretching along its own direction only (NOTES.md section 4, issue 1).

The test: a straight strip of the wall's thickness, meshed exactly like the
fin (edge points every h = ds/sqrt(mult), Halton fill, the same boundary gap,
the same spring formula, the same target clamp), clamped at x = 0 with target
points, loaded by a force F spread evenly over its free end. The static
deflection at scale 1 is compared with Timoshenko beam theory (bending plus
shear):

    tip deflection = F L^3 / (3 E' I) + F L / (kappa G A)
    I = t^3 / 12, A = t (per unit depth), kappa = 5/6, G = E / (2 (1 + nu))
    E' = E (plane stress) or E / (1 - nu^2) (plane strain)

Deflection is exactly proportional to 1 / scale (linear solve), so
scale = mesh deflection at scale 1 / reference deflection. Two estimates:
  scale_tip  from the tip deflection (mean over the end face)
  scale_fit  least-squares fit of the centreline deflection over
             FIT_FROM*L <= s <= L (leaves out the region next to the clamp,
             where the clamp's local deformation is not beam-like)

Optional IB2d cross-check (--ib2d-check): IB2d's own spring + target force
routines at the solved shape must balance the load; the residual must fall
in proportion to the load (only the geometric nonlinearity a linear solve
leaves out remains), as on the wing.

The static problem is the linear stiffness matrix of the springs (k n n^T
each) and target tethers, solved sparse. Units cancel: IB2d stiffnesses and
force densities are both the physical values / ds, so only h matters, not ds.

Defaults (decided 2026-10-18, EXPERIMENT_PLAN.md decisions log): two strips,
the fin's long side wall (138.4 mm, L/t = 35: almost pure bending) and its
shortest ray (16.9 mm, L/t = 4.3: bending plus ~6% shear), both centreline
lengths; plane strain (the wall is 27 mm deep for 4 mm thick, so it bends
like a plate, and the 2D model is per unit depth); nu = 0.499 (silicone).

Run (light: seconds; still run it on a compute node, see CLUSTER.md):
    python verify_structure.py [--mults 1,4,...] [--seeds none,0,1] [options]
Writes one CSV row per (length, multiplier, seed) to --out.
"""

import argparse
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "common"))
import mesh_utils as mu
from ib_units import ib_ds

sys.path.insert(0, str(HERE))
import FinRay_Geom as FR

KAPPA = 5.0 / 6.0       # Timoshenko shear coefficient, rectangular section
FIT_FROM = 0.1          # centreline fit over FIT_FROM*L <= s <= L
TIP_DEFLECTION = 1e-3   # load set so the bending tip deflection is this fraction of L (linear)


# --------------------------------------------------------------------------
# Linear static solve (as cases/NACA0020_Wing/verify_structure.py)
# --------------------------------------------------------------------------

def stiffness_matrix(x, y, connections, edge_k, target_idx, k_target):
    '''Sparse linear stiffness (IB2d units) of the springs (k n n^T each) and
    the target tethers, at the rest configuration. DOF order x0, y0, x1, ...'''
    ij = np.array([(i, j) for i, j, _ in connections])
    k = np.array([edge_k[(i, j)] for i, j, _ in connections])
    d = np.column_stack([x[ij[:, 1]] - x[ij[:, 0]], y[ij[:, 1]] - y[ij[:, 0]]])
    n = d / np.linalg.norm(d, axis=1)[:, None]
    rows, cols, vals = [], [], []
    for a in range(2):
        for b in range(2):
            kab = k * n[:, a] * n[:, b]
            for p, q, s in ((0, 0, 1), (1, 1, 1), (0, 1, -1), (1, 0, -1)):
                rows.append(2 * ij[:, p] + a); cols.append(2 * ij[:, q] + b); vals.append(s * kab)
    t = np.asarray(target_idx)
    for a in range(2):
        rows.append(2 * t + a); cols.append(2 * t + a); vals.append(np.full(len(t), k_target))
    N = 2 * len(x)
    return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                         shape=(N, N))


def ib2d_residual(x, y, connections, edge_k, target_idx, k_target, ds, ux, uy, fy, eps):
    '''max |IB2d spring + target forces + load| / max |load| at the solved
    shape, with displacement and load scaled by eps.'''
    sys.path.insert(0, str(HERE.parents[1] / "external" / "IB2d" / "pyIB2d" / "IBM_Blackbox"))
    from please_Find_Lagrangian_Forces_On_Eulerian_grid import (
        give_Me_Spring_Lagrangian_Force_Densities as spring_f,
        give_Me_Target_Lagrangian_Force_Densities as target_f)
    springs = np.array([(i, j, edge_k[(i, j)], np.hypot(x[j] - x[i], y[j] - y[i]), 1.0)
                        for i, j, _ in connections])
    t = np.asarray(target_idx)
    targets = np.column_stack([t, x[t], y[t], np.full(len(t), k_target)])
    xs, ys = x + eps * ux, y + eps * uy
    big = 1e3   # domain size for IB2d's periodic wrap: large = no wrap
    fx1, fy1 = spring_f(ds, len(xs), xs, ys, springs, big, big)
    fx2, fy2 = target_f(ds, xs, ys, targets, big, big)
    r = np.concatenate([fx1 + fx2, fy1 + fy2 + eps * fy])
    return float(np.abs(r).max() / np.abs(eps * fy).max())


# --------------------------------------------------------------------------
# Strip test
# --------------------------------------------------------------------------

def fin_boundary_gap(wall_thickness):
    '''The fill's min_boundary_gap in the real fin (fill_polygon_with_halton's
    default for the whole tail polygon), so the strip is meshed by the same rule.'''
    outer, windows = FR._tail_outline(wall_thickness)
    loops = [outer] + windows
    area = mu._shoelace_area(outer) - sum(mu._shoelace_area(w) for w in windows)
    perimeter = sum(np.linalg.norm(np.roll(l, -1, axis=0) - l, axis=1).sum() for l in loops)
    return 0.25 * area / perimeter


def timoshenko(s, L, F, E, t, nu, plane):
    '''Deflection of a clamped beam under an end force F at stations s, split
    into (bending, shear).'''
    Ep = E / (1 - nu ** 2) if plane == "strain" else E
    EI = Ep * t ** 3 / 12.0
    G = E / (2 * (1 + nu))
    return F / EI * (L * s ** 2 / 2 - s ** 3 / 6), F * s / (KAPPA * G * t)


def strip_test(L, t, ds, mult, seed, E, nu, plane, target_factor, gap, ib2d_check=False):
    rect = np.array([[0.0, 0.0], [L, 0.0], [L, t], [0.0, t]])
    h = ds / np.sqrt(max(mult, 1.0))
    loop = mu.resample_loop(rect, h)
    fill = mu.fill_polygon_with_halton(rect, ds, density_multiplier=mult,
                                       min_boundary_gap=gap, seed=seed)
    pts, _, tris = mu.triangulate_with_holes(loop, [], loop, fill)
    x, y = pts[:, 0], pts[:, 1]
    conn, ek = mu.tributary_springs(x, y, tris, E, ds, 1.0)
    clamp = np.where(x < 1e-12)[0]
    kT = target_factor * mu.max_incident_stiffness(clamp.tolist(), ek)

    # end force F (per unit depth) spread evenly over the end face, as IB2d force densities
    Ep = E / (1 - nu ** 2) if plane == "strain" else E
    F = TIP_DEFLECTION * L * 3 * Ep * t ** 3 / 12.0 / L ** 3
    end = np.where(x > L - 1e-12)[0]
    ye = np.sort(y[end]); order = np.argsort(y[end])
    edges = np.concatenate([[0.0], 0.5 * (ye[1:] + ye[:-1]), [t]])
    trib = np.empty(len(end)); trib[order] = np.diff(edges)
    fy = np.zeros(len(x)); fy[end] = F * trib / t / ds

    K = stiffness_matrix(x, y, conn, ek, clamp, kT)
    f = np.zeros(2 * len(x)); f[1::2] = fy
    u = spla.spsolve(K.tocsc(), f)
    ux, uy = u[0::2], u[1::2]

    import matplotlib.tri as mtri
    s = np.linspace(FIT_FROM, 1.0, 19) * L
    s[-1] = L * (1 - 1e-9)   # stay inside the triangulation
    u_cl = np.asarray(mtri.LinearTriInterpolator(mtri.Triangulation(x, y, tris), uy)(s, np.full(len(s), t / 2)))
    bend, shear = timoshenko(s, L, F, E, t, nu, plane)
    ref = bend + shear
    tip_mesh = float(np.sum(uy[end] * trib) / t)
    tb, ts = timoshenko(np.array([L]), L, F, E, t, nu, plane)
    tip_ref = float(tb[0] + ts[0])

    row = dict(mult=mult, seed="none" if seed is None else seed, h_mm=1e3 * h,
               points_across=t / h, n_points=len(x), n_springs=len(conn),
               scale_tip=tip_mesh / tip_ref,
               scale_fit=float(np.sum(u_cl * u_cl) / np.sum(u_cl * ref)),
               shear_share=float(ts[0] / tip_ref))
    if ib2d_check:
        row["ib2d_residual"] = ";".join(
            f"{ib2d_residual(x, y, conn, ek, clamp, kT, ds, ux, uy, fy, e):.1e}" for e in (1.0, 1e-2, 1e-4))
    return row


def main():
    params = json.loads((HERE / "case_params.json").read_text())
    g, grid = params["geom"], params["input2d"]
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--lengths", default="0.1384,0.0169",
                    help="strip lengths L [m], comma-separated (default: side wall, shortest ray)")
    ap.add_argument("--plane", choices=("stress", "strain"), default="strain",
                    help="reference: plane stress (E' = E) or plane strain (E' = E/(1-nu^2), default)")
    ap.add_argument("--thickness", type=float, default=g["wall_thickness"],
                    help="strip thickness t [m] (default: geom.wall_thickness)")
    ap.add_argument("--nu", type=float, default=0.499, help="material Poisson's ratio for the reference")
    ap.add_argument("--nx", type=int, default=int(grid["Nx"]), help="grid Nx, sets ds (default: input2d.Nx)")
    ap.add_argument("--mults", default="1,4,16,64,256,1024", help="density multipliers, comma-separated")
    ap.add_argument("--seeds", default="none",
                    help="Halton seeds, comma-separated; 'none' = the unscrambled fill used by default")
    ap.add_argument("--ib2d-check", action="store_true", help="IB2d force cross-check on every row (slower)")
    ap.add_argument("--out", default="verify_structure.csv", help="CSV output path")
    a = ap.parse_args()

    E, t = g["E_material"], a.thickness
    ds = ib_ds(float(grid["Lx"]), a.nx)
    gap = fin_boundary_gap(g["wall_thickness"])
    mults = [float(m) for m in a.mults.split(",")]
    seeds = [None if s.strip().lower() == "none" else int(s) for s in a.seeds.split(",")]
    print(f"t = {1e3 * t:.2f} mm, Nx = {a.nx} (ds = {1e3 * ds:.3f} mm), plane {a.plane}, "
          f"nu = {a.nu}, E = {E:.3g} Pa, boundary gap {1e3 * gap:.3f} mm")
    rows = []
    for L in (float(v) for v in a.lengths.split(",")):
        print(f"\nStrip L = {1e3 * L:.1f} mm (L/t = {L / t:.1f})")
        print("   mult   seed   h[mm]  pts/t   points  scale_tip  scale_fit  shear%  time[s]")
        for mult in mults:
            for seed in seeds:
                t0 = time.perf_counter()
                r = strip_test(L, t, ds, mult, seed, E, a.nu, a.plane, g["target_stiffness_factor"],
                               gap, a.ib2d_check)
                r.update(length=L, thickness=t, nx=a.nx, plane=a.plane, nu=a.nu)
                rows.append(r)
                print(f"   {mult:5g}  {r['seed']:>5}  {r['h_mm']:6.3f}  {r['points_across']:5.1f}  "
                      f"{r['n_points']:7d}  {r['scale_tip']:9.4f}  {r['scale_fit']:9.4f}  "
                      f"{100 * r['shear_share']:5.1f}  {time.perf_counter() - t0:7.1f}"
                      + (f"   IB2d residual {r['ib2d_residual']}" if a.ib2d_check else ""))
        if len(seeds) > 1:
            print("   mult   scale_tip mean +- std   scale_fit mean +- std   (over seeds)")
            for mult in mults:
                st = np.array([r["scale_tip"] for r in rows if r["mult"] == mult and r["length"] == L])
                sf = np.array([r["scale_fit"] for r in rows if r["mult"] == mult and r["length"] == L])
                print(f"   {mult:5g}   {st.mean():.4f} +- {st.std(ddof=1):.4f}        "
                      f"{sf.mean():.4f} +- {sf.std(ddof=1):.4f}")
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
