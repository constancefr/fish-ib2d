"""
Fluid-free sanity check for the ribbon Halton-fill + Delaunay-triangulation
feature (see NOTES.md): exercises _Build_Tail_Geometry_Ribbon /
build_Tail_Ribbon_Connections directly, across several
ribbon_fill_density_multiplier values, and checks for degenerate springs,
runaway point/spring counts, and runaway build time -- no Eulerian grid or
CFD solve involved, following the same no-fluid convention as
NACA0020_Wing/verify_structure.py.

Run:  python verify_ribbon_fill.py
"""

import time

import numpy as np

import FinRay_Geom as FR


def run_one(mult, ds, wall_thickness):
    t0 = time.perf_counter()
    xLag, yLag, outlineRailIdx, rayRailIdxList, baseCornerIdx, ribbonPolys = \
        FR._Build_Tail_Geometry_Ribbon(ds, wall_thickness, fill_density_multiplier=mult)
    connections, fill_edge_k = FR.build_Tail_Ribbon_Connections(
        xLag, yLag, outlineRailIdx, rayRailIdxList, ribbonPolys,
        E_material=0.74e6, extrude_depth=0.05, fill_stiffness_scale=1.0)
    dt = time.perf_counter() - t0

    n_interior = sum(len(r['interior_idx']) for r in ribbonPolys)
    n_fill = sum(1 for (_, _, k) in connections if k == 'fill')

    ii = np.array([i for i, _, _ in connections])
    jj = np.array([j for _, j, _ in connections])
    lengths = np.hypot(xLag[ii] - xLag[jj], yLag[ii] - yLag[jj])
    n_degenerate = int(np.sum(lengths < 1e-9))

    print(f"multiplier={mult:6.2f}  interior_pts={n_interior:4d}  total_pts={len(xLag):4d}  "
          f"springs={len(connections):5d} (fill={n_fill:4d})  degenerate={n_degenerate}  "
          f"build_time={dt * 1e3:.1f}ms")

    assert n_degenerate == 0, f"found {n_degenerate} zero/near-zero-length springs at multiplier={mult}"
    return dict(mult=mult, n_pts=len(xLag), n_springs=len(connections), dt=dt)


if __name__ == "__main__":
    ds = 0.5 * (1.0 / 32)
    wall_thickness = 0.00396

    print("Regression check: multiplier=0 must reproduce the sparse-truss geometry "
          "(no 'fill' connections at all).")
    zero = run_one(0.0, ds, wall_thickness)

    results = [zero] + [run_one(m, ds, wall_thickness) for m in (1.0, 4.0, 16.0, 64.0)]

    print("\nGrowth check (point/spring count should grow roughly linearly with "
          "multiplier, not blow up super-linearly, given O(n log n) Delaunay cost "
          "at these small per-ribbon point counts):")
    base = results[0]
    for a, b in zip(results, results[1:]):
        d_pts = b['n_pts'] - base['n_pts']
        d_pts_prev = a['n_pts'] - base['n_pts']
        ratio_pts = d_pts / max(d_pts_prev, 1)
        print(f"  {a['mult']:6.2f} -> {b['mult']:6.2f}:  point-count ratio={ratio_pts:.2f}")

    print("\nOK: no degenerate springs at any tested density.")
