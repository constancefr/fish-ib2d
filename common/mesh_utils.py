"""
Shape-agnostic Halton-fill + constrained-Delaunay utilities for turning a
polygon (with optional holes) into a "solid" 2D Lagrangian point mesh:
interior points + a triangulation-derived spring edge list.

Nothing here knows about ribbons, rays, tails, or wings -- callers hand in a
polygon boundary (and, for an annular/hollow region, hole boundaries) and get
back interior points, mesh connectivity and tributary-area springs. Used by
cases/FinRay/FinRay_Geom.py and cases/NACA0020_Wing/Wing_Geom.py.
"""

import numpy as np
from matplotlib.path import Path as _MplPath
import meshpy.triangle as triangle
from scipy.stats import qmc

from ib_units import to_ib


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

def _dist_to_segments(p, segments):
    '''Distance from each point in p (n,2) to the nearest of the (a, b) segments.'''
    d = np.full(len(p), np.inf)
    for a, b in segments:
        t = np.clip((p - a) @ (b - a) / max((b - a) @ (b - a), 1e-30), 0.0, 1.0)
        d = np.minimum(d, np.linalg.norm(p - (a + t[:, None] * (b - a)), axis=1))
    return d


def _loop_segments(loop):
    '''The (a, b) sides of a closed loop (first vertex NOT repeated).'''
    return list(zip(loop, np.roll(loop, -1, axis=0)))


def fill_polygon_with_halton(exterior, spacing, inner_loops=None,
                              density_multiplier=1.0, min_boundary_gap=None,
                              keep_out_segments=None):
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
    keep_out_segments: extra (a, b) segments inside the polygon that fill
        points also keep min_boundary_gap from, e.g. an internal line that
        will be passed to triangulate_with_holes as extra_segments.

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

    segments = [s for l in loops for s in _loop_segments(l)]
    segments += [(np.asarray(a, dtype=float), np.asarray(b, dtype=float))
                 for a, b in (keep_out_segments or [])]
    path_ext = closed_path(exterior)
    path_holes = [closed_path(h) for h in inner_loops]
    pts = np.empty((0, 2))
    for _ in range(8):  # safety cap; one round is normally enough
        cand = qmc.scale(sampler.random(batch), lo, hi)  # continues the sequence each round
        inside = path_ext.contains_points(cand)
        for ph in path_holes:
            inside &= ~ph.contains_points(cand)
        cand = cand[inside]
        cand = cand[_dist_to_segments(cand, segments) >= min_boundary_gap]
        pts = np.vstack([pts, cand])
        if len(pts) >= n:
            return pts[:n]

    print(f"[mesh_utils] fill_polygon_with_halton: placed {len(pts)}/{n} points "
          f"(polygon too thin for this spacing/min_boundary_gap); continuing with fewer.")
    return pts


# --- Constrained Delaunay triangulation -> edge list -------------------------

def triangulate_with_holes(exterior, inner_loops, boundary_pts, interior_pts, extra_segments=None):
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
    extra_segments: further (i, j) local-index pairs kept as triangle edges
        (e.g. an internal line no triangle may cross). Points on such a line
        that are not loop points go in interior_pts.

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
    segments += [(int(i), int(j)) for i, j in (extra_segments or [])]

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


# --- Springs from a triangulation --------------------------------------------

def tributary_springs(xLag, yLag, triangles, E, ds, stiffness_scale=1.0):
    '''
    One spring per edge of the triangulation, each with the same
    tributary-area stiffness:

        k = to_ib(stiffness_scale * E * A_tributary / L**2, ds)
          = stiffness_scale * E * A_tributary / (L**2 * ds)

    i.e. a physical stiffness per unit depth [N/m per m], divided by IB2d's
    force-spreading length ds (see common/ib_units.py).

    A_tributary is 1/3 of the area of each triangle sharing that edge, so
    the springs hold every triangle's area exactly once, at any sampling
    density. This does NOT make the network behave like a sheet of modulus
    E: a random Halton/Delaunay network with stiffness_scale = 1 has an
    effective Young's modulus of only ~0.33-0.36 E (exactly E/3 for an
    equilateral mesh) and a Poisson's ratio of ~0.3, because each spring
    resists stretching along its own direction only. stiffness_scale (the
    single calibration multiplier) therefore has to be measured, e.g. by a
    static bending test (see cases/NACA0020_Wing/verify_structure.py).
    An outline edge borders one triangle only, so it gets roughly half the
    stiffness of a similar interior edge.

    kind (plotting only): 'boundary' for outline edges (one triangle),
    'fill' for every other edge.

    Returns (connections, edge_k): connections is a list of (i, j, kind);
    edge_k maps each connection's (min(i,j), max(i,j)) pair to its k.
    '''
    edge_area, edge_ntri = {}, {}
    for tri in triangles:
        x, y = xLag[tri], yLag[tri]
        area = 0.5 * abs((x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]))
        for a, b in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            pair = (int(min(a, b)), int(max(a, b)))
            edge_area[pair] = edge_area.get(pair, 0.0) + area / 3.0
            edge_ntri[pair] = edge_ntri.get(pair, 0) + 1

    connections = [(i, j, 'boundary' if edge_ntri[(i, j)] == 1 else 'fill')
                   for (i, j) in sorted(edge_area)]
    edge_k = {}
    for (i, j), A in edge_area.items():
        L = np.hypot(xLag[i] - xLag[j], yLag[i] - yLag[j])
        edge_k[(i, j)] = to_ib(stiffness_scale * E * A / L ** 2, ds)
    return connections, edge_k


def max_incident_stiffness(point_indices, edge_k):
    '''
    Largest total spring stiffness attached to any one of point_indices,
    i.e. sum of k over every spring touching that point, maxed over points.
    Used as the reference magnitude for k_Target (target points must
    overpower the springs they pull on, but not by so much that the explicit
    time stepping blows up).

    edge_k: dict (min(i,j), max(i,j)) -> k, from tributary_springs.
    '''
    incident = {p: 0.0 for p in point_indices}
    for (i, j), k in edge_k.items():
        if i in incident:
            incident[i] += k
        if j in incident:
            incident[j] += k
    return max(incident.values())
