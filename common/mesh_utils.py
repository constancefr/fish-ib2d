"""
Shape-agnostic Halton-fill + constrained-Delaunay utilities for turning a
polygon (with optional holes) into a "solid" 2D Lagrangian point mesh:
interior points + a triangulation-derived spring edge list.

Nothing here knows about ribbons, rays, tails, or wings -- callers hand in a
polygon boundary (and, for an annular/hollow region, hole boundaries) and get
back interior points plus mesh connectivity. Reusable by any
cases/*/*_Geom.py builder (currently FinRay_Geom.py; potentially
NACA0020_Wing/Wing_Geom.py or a fully-filled-shape example later).
"""

import numpy as np
from matplotlib.path import Path as _MplPath
import meshpy.triangle as triangle
from scipy.stats import qmc


def _shoelace_area(loop):
    '''
    Area of a simple polygon given as its vertex list (first vertex NOT repeated at the end).
    Uses the shoelace formula.
    '''
    loop = np.asarray(loop, dtype=float)
    x, y = loop[:, 0], loop[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def closed_path(loop):
    '''
    matplotlib Path for a closed polygon given as its vertex list (first
    vertex NOT repeated at the end). Path(loop, closed=True) on its own
    treats the LAST vertex as a CLOSEPOLY placeholder and ignores its
    coordinates, silently dropping that corner from every contains_point
    test -- so the first vertex is appended here to stand in for it.
    '''
    loop = np.asarray(loop, dtype=float)
    return _MplPath(np.vstack([loop, loop[:1]]), closed=True)


def resample_loop(loop, spacing):
    '''
    Points every ~spacing along a closed polygon (first vertex NOT repeated),
    keeping its vertices: each side gets max(1, round(length / spacing))
    equal segments.
    '''
    loop = np.asarray(loop, dtype=float)
    out = []
    for a, b in zip(loop, np.roll(loop, -1, axis=0)):
        n = max(1, round(np.linalg.norm(b - a) / spacing))
        out += [a + (b - a) * k / n for k in range(n)]
    return np.array(out)


# --- Interior fill sampler ---------------------------------------------------

def _dist_to_loops(p, loops):
    '''Distance from each point in p (n,2) to the nearest segment of any closed loop.'''
    d = np.full(len(p), np.inf)
    for loop in loops:
        for a, b in zip(loop, np.roll(loop, -1, axis=0)):
            t = np.clip((p - a) @ (b - a) / max((b - a) @ (b - a), 1e-30), 0.0, 1.0)
            d = np.minimum(d, np.linalg.norm(p - (a + t[:, None] * (b - a)), axis=1))
    return d


def fill_polygon_with_halton(exterior, spacing, inner_loops=None,
                              density_multiplier=1.0, min_boundary_gap=None):
    '''
    Fills the polygon (exterior minus inner_loops) with Halton-sampled
    interior points: round(density_multiplier * area / spacing**2) of them,
    i.e. on average density_multiplier points per spacing x spacing cell.

    exterior: (n, 2) simple closed polygon (do not repeat the first vertex).
    inner_loops: list of (m_k, 2) closed polygons (holes) inside `exterior`.
    density_multiplier: 0 disables fill (returns no points).
    min_boundary_gap: no point closer than this to any edge (exterior or
        hole), so fill points don't sit almost on a boundary segment and
        create near-flat triangles. Defaults to 0.25 * area / perimeter,
        about 1/8 of a thin strip's width (not tied to `spacing`: for a strip
        narrower than spacing, a spacing-based gap would reject everything).

    Candidates are drawn in the polygon's x/y bounding box and kept only if
    inside and far enough from the boundary (rejection sampling). Keeping
    the first n in draw order still covers the polygon evenly. Prints a
    warning and returns fewer points if the polygon is too thin to fit n.
    '''
    exterior = np.asarray(exterior, dtype=float)
    inner_loops = [np.asarray(h, dtype=float) for h in (inner_loops or [])]
    loops = [exterior] + inner_loops

    area = _shoelace_area(exterior) - sum(_shoelace_area(h) for h in inner_loops)
    n = int(round(density_multiplier * area / spacing ** 2))
    if n <= 0:
        return np.empty((0, 2))
    if min_boundary_gap is None:
        perimeter = sum(np.linalg.norm(np.roll(l, -1, axis=0) - l, axis=1).sum() for l in loops)
        min_boundary_gap = 0.25 * area / perimeter

    lo, hi = exterior.min(axis=0), exterior.max(axis=0)
    batch = int(np.ceil(2 * n * np.prod(hi - lo) / area))  # ~2x the candidates we expect to need
    sampler = qmc.Halton(d=2, scramble=False)
    sampler.fast_forward(1)  # index 0 is the origin (a degenerate point)

    path_ext = closed_path(exterior)
    path_holes = [closed_path(h) for h in inner_loops]
    pts = np.empty((0, 2))
    for _ in range(8):  # safety cap; one round is normally enough
        cand = qmc.scale(sampler.random(batch), lo, hi)  # continues the sequence each round
        inside = path_ext.contains_points(cand)
        for ph in path_holes:
            inside &= ~ph.contains_points(cand)
        cand = cand[inside]
        cand = cand[_dist_to_loops(cand, loops) >= min_boundary_gap]
        pts = np.vstack([pts, cand])
        if len(pts) >= n:
            return pts[:n]

    print(f"[mesh_utils] fill_polygon_with_halton: placed {len(pts)}/{n} points "
          f"(polygon too thin for this spacing/min_boundary_gap); continuing with fewer.")
    return pts


# --- Constrained Delaunay triangulation -> edge list -------------------------

def triangulate_with_holes(exterior, inner_loops, boundary_pts, interior_pts):
    '''
    Constrained Delaunay triangulation of boundary_pts ++ interior_pts with
    Triangle (via meshpy): every `exterior`/`inner_loops` segment is kept as
    a triangle edge, and triangles outside `exterior` or inside a hole are
    dropped. No points are added.

    boundary_pts must be exactly the points of `exterior` followed by the
    points of each `inner_loops[k]`, in that same order -- callers use local
    indices 0..nb-1 to map back to their own global point indices.
    interior_pts occupy local indices nb..nb+ni-1. Each hole is marked by
    its vertex mean, so holes must contain it (true for convex holes).

    Returns (all_pts, edges, triangles):
      all_pts: (nb+ni, 2) = vstack(boundary_pts, interior_pts)
      edges: sorted list of unique (i, j) local-index pairs, i < j
      triangles: (nt, 3) array of local vertex indices
    '''
    inner_loops = [np.asarray(h, dtype=float) for h in (inner_loops or [])]
    interior_pts = np.asarray(interior_pts, dtype=float).reshape(-1, 2)
    all_pts = np.vstack([np.asarray(boundary_pts, dtype=float), interior_pts])

    segments, start = [], 0
    for n in [len(exterior)] + [len(h) for h in inner_loops]:
        segments += [(start + k, start + (k + 1) % n) for k in range(n)]
        start += n

    hole_marks = [h.mean(axis=0) for h in inner_loops]
    assert all(closed_path(h).contains_point(m) for h, m in zip(inner_loops, hole_marks)), \
        "a hole does not contain its vertex mean; give it an explicit marker point"

    info = triangle.MeshInfo()
    info.set_points(all_pts)
    info.set_facets(segments)
    info.set_holes(hole_marks)
    mesh = triangle.build(info, quality_meshing=False,  # no added points
                          allow_boundary_steiner=False, allow_volume_steiner=False)
    assert len(mesh.points) == len(all_pts), "Triangle added or merged points"

    triangles = np.array(mesh.elements, dtype=int).reshape(-1, 3)
    edges = sorted({(min(a, b), max(a, b))
                    for t in triangles for a, b in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0]))})
    return all_pts, edges, triangles
