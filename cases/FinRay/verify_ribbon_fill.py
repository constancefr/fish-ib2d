"""
Fluid-free sanity check for the tail mesh (one outline + Halton fill + one
constrained triangulation, see NOTES.md): exercises build_Tail_Mesh /
build_Tail_Springs directly, across several
ribbon_fill_density_multiplier values, and checks for degenerate springs,
the material budget, runaway point/spring counts, and runaway build time --
no Eulerian grid or CFD solve involved, following the same no-fluid
convention as NACA0020_Wing/verify_structure.py.

Run:  python verify_ribbon_fill.py
"""

import time

import numpy as np

import FinRay_Geom as FR

E = 0.74e6


def run_one(mult, ds, wall_thickness):
    t0 = time.perf_counter()
    xLag, yLag, triangles, n_boundary, _, outer, windows = \
        FR.build_Tail_Mesh(ds, wall_thickness, ds, mult)
    connections, edge_k = FR.build_Tail_Springs(
        xLag, yLag, triangles, E_material=E, ds=ds, stiffness_scale=1.0)
    dt = time.perf_counter() - t0

    n_interior = len(xLag) - n_boundary
    n_fill = sum(1 for (_, _, k) in connections if k == 'fill')

    ii = np.array([i for i, _, _ in connections])
    jj = np.array([j for _, j, _ in connections])
    lengths = np.hypot(xLag[ii] - xLag[jj], yLag[ii] - yLag[jj])
    n_degenerate = int(np.sum(lengths < 1e-9))

    print(f"multiplier={mult:6.2f}  edge_pts={n_boundary:4d}  interior_pts={n_interior:4d}  "
          f"total_pts={len(xLag):4d}  springs={len(connections):5d} (fill={n_fill:4d})  "
          f"degenerate={n_degenerate}  build_time={dt * 1e3:.1f}ms")

    assert n_degenerate == 0, f"found {n_degenerate} zero/near-zero-length springs at multiplier={mult}"

    # Material budget in IB2d units: k = E*A_trib/(L^2*ds), so the springs'
    # sum of k*ds*L^2/E must equal the tail's area (material counted once).
    k = np.array([edge_k[(min(i, j), max(i, j))] for i, j, _ in connections])
    budget = np.sum(k * ds * lengths**2 / E)
    area = FR.mu._shoelace_area(outer) - sum(FR.mu._shoelace_area(w) for w in windows)
    print(f"    material budget = {100 * budget / area:.4f}% of tail area ({area * 1e6:.1f} mm^2)")
    # A unit error would be off by orders of magnitude; anything else means
    # the triangles overlap or leave gaps.
    assert abs(budget / area - 1) < 1e-2, f"material budget off at multiplier={mult}"
    if abs(budget / area - 1) > 1e-9:
        print("    WARNING: triangles do not tile the tail exactly (overlap or gap)")
    return dict(mult=mult, n_pts=len(xLag), n_springs=len(connections), dt=dt)


if __name__ == "__main__":
    ds = 0.5 * (1.0 / 32)
    wall_thickness = 0.00396

    print("multiplier=0: edge points only (every ~ds), no interior points ('fill' springs are then "
          "just the cross-wall triangulation edges).")
    zero = run_one(0.0, ds, wall_thickness)

    results = [zero] + [run_one(m, ds, wall_thickness) for m in (1.0, 4.0, 16.0, 64.0)]

    print("\nGrowth check (point/spring count should grow roughly linearly with "
          "multiplier, not blow up super-linearly, given O(n log n) Delaunay cost "
          "at these point counts):")
    base = results[0]
    for a, b in zip(results, results[1:]):
        d_pts = b['n_pts'] - base['n_pts']
        d_pts_prev = a['n_pts'] - base['n_pts']
        ratio_pts = d_pts / max(d_pts_prev, 1)
        print(f"  {a['mult']:6.2f} -> {b['mult']:6.2f}:  point-count ratio={ratio_pts:.2f}")

    print("\nOK: no degenerate springs and material budget within 1% at every tested density.")
