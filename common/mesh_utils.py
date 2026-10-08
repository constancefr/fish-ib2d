"""
Shape-agnostic Halton-fill + Delaunay-triangulation utilities for turning a
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

try:
    from scipy.spatial import Delaunay as _Delaunay
    _SCIPY_IMPORT_ERROR = None
except ImportError as _e:
    _Delaunay = None
    _SCIPY_IMPORT_ERROR = _e


# --- Halton sequence --------------------------------------------------------

def _van_der_corput(idxs, base):
    '''Radical-inverse of each integer in idxs, in the given base.'''
    result = np.zeros(idxs.shape, dtype=float)
    i = idxs.copy()
    f = 1.0 / base
    while np.any(i > 0):
        result += f * (i % base)
        i //= base
        f /= base
    return result


def halton(n, dims=2, bases=(2, 3), start=1):
    '''
    Van der Corput / Halton low-discrepancy sequence.

    n: number of points to generate
    dims: dimensionality (<= len(bases))
    bases: one coprime base per dimension (default (2, 3): standard 2D Halton)
    start: 1-based starting index. Index 0 maps to the origin in every base
        (a degenerate point), so the default start=1 skips it.

    Returns an (n, dims) array of points in [0, 1)^dims.
    '''
    idxs = np.arange(start, start + n, dtype=np.int64)
    pts = np.empty((n, dims), dtype=float)
    for d in range(dims):
        pts[:, d] = _van_der_corput(idxs, bases[d])
    return pts


# --- Oriented bounding box ---------------------------------------------------

def oriented_bbox(vertices):
    '''
    PCA-based oriented bounding box of a point set (e.g. a polygon's
    vertices). Used to draw efficient Halton candidates for elongated /
    rotated shapes, where a plain axis-aligned bbox would waste most of its
    candidates on empty space.

    Returns (center, u, v, half_u, half_v):
      u, v: orthonormal axes (u = direction of max vertex-position variance)
      half_u, half_v: half-extents of the vertex projections along u, v

    A point p is inside the box iff
      abs(dot(p - center, u)) <= half_u and abs(dot(p - center, v)) <= half_v
    '''
    pts = np.asarray(vertices, dtype=float)
    centroid = pts.mean(axis=0)
    centered = pts - centroid
    cov = centered.T @ centered / len(pts)
    eigvals, eigvecs = np.linalg.eigh(cov)  # ascending order
    order = np.argsort(eigvals)[::-1]
    u = eigvecs[:, order[0]]
    v = eigvecs[:, order[1]]

    proj_u = centered @ u
    proj_v = centered @ v
    center = centroid + 0.5 * (proj_u.min() + proj_u.max()) * u \
                       + 0.5 * (proj_v.min() + proj_v.max()) * v
    half_u = 0.5 * (proj_u.max() - proj_u.min())
    half_v = 0.5 * (proj_v.max() - proj_v.min())
    return center, u, v, half_u, half_v


def _shoelace_area(loop):
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


def _perimeter(loop):
    loop = np.asarray(loop, dtype=float)
    return float(np.sum(np.linalg.norm(np.roll(loop, -1, axis=0) - loop, axis=1)))


# --- Point-to-polygon-boundary distance (for min_boundary_gap) -------------

def _point_segment_dist(p, a, b):
    '''Distance from point(s) p (n,2) to the segment a->b.'''
    ab = b - a
    denom = np.dot(ab, ab)
    if denom < 1e-30:
        return np.linalg.norm(p - a, axis=1)
    t = np.clip(((p - a) @ ab) / denom, 0.0, 1.0)
    proj = a + t[:, None] * ab
    return np.linalg.norm(p - proj, axis=1)


def _min_dist_to_loop(p, loop):
    loop = np.asarray(loop, dtype=float)
    n = len(loop)
    dists = np.full(len(p), np.inf)
    for k in range(n):
        a, b = loop[k], loop[(k + 1) % n]
        dists = np.minimum(dists, _point_segment_dist(p, a, b))
    return dists


# --- Point-to-polyline projection (for welding one ribbon to another) -------

def project_point_to_polyline(p, polyline, closed=True):
    '''
    The closest point to p that lies ON the polyline itself (not just on one
    of its vertices) -- e.g. for finding where a ray ribbon's edge actually
    meets the outline ribbon's boundary, instead of snapping to whichever
    existing discretized boundary point happens to be nearest.

    p: (2,) query point
    polyline: (n,2) ordered vertices
    closed: True if the last vertex implicitly connects back to the first

    Returns (foot_xy, seg_idx, t):
      foot_xy: the closest point on the polyline
      seg_idx: index i such that foot_xy lies on the segment
          polyline[i] -> polyline[(i+1) % n]
      t: how far along that segment, in [0, 1] (0 = polyline[seg_idx],
          1 = polyline[(seg_idx+1) % n])
    '''
    p = np.asarray(p, dtype=float)
    polyline = np.asarray(polyline, dtype=float)
    n = len(polyline)
    stop = n if closed else n - 1

    best_d2 = np.inf
    best = None
    for i in range(stop):
        a = polyline[i]
        b = polyline[(i + 1) % n]
        ab = b - a
        denom = np.dot(ab, ab)
        t = 0.0 if denom < 1e-30 else float(np.clip(np.dot(p - a, ab) / denom, 0.0, 1.0))
        foot = a + t * ab
        d2 = float(np.dot(p - foot, p - foot))
        if d2 < best_d2:
            best_d2 = d2
            best = (foot, i, t)
    return best


def _far_from_boundary(p, exterior, inner_loops, min_gap):
    dists = _min_dist_to_loop(p, exterior)
    for h in inner_loops:
        dists = np.minimum(dists, _min_dist_to_loop(p, h))
    return dists >= min_gap


# --- Interior fill sampler ---------------------------------------------------

def fill_polygon_with_halton(exterior, inner_loops=None, spacing=None,
                              n_target=None, density_multiplier=1.0,
                              min_boundary_gap=None, max_rounds=8,
                              bases=(2, 3)):
    '''
    Fills the polygon (exterior minus inner_loops) with Halton-sampled
    interior points at area density ~ density_multiplier / spacing**2 (on
    average one point per spacing x spacing cell), or an explicit n_target
    count.

    exterior: (n, 2) simple closed polygon (do not repeat the first vertex).
    inner_loops: list of (m_k, 2) simple closed polygons (holes), each
        entirely inside `exterior` and mutually non-overlapping. None/[]
        for a simply-connected polygon.
    spacing: target point spacing. Ignored if n_target is given.
    n_target: overrides the density-derived count with an exact number.
    density_multiplier: scales the density-derived point count. 0 disables
        fill entirely (returns no points).
    min_boundary_gap: minimum allowed distance from any accepted point to
        the polygon boundary (exterior or any hole edge). Defaults to
        0.25 * (poly_area / perimeter), an area/perimeter estimate of the
        polygon's local half-width that stays meaningful for a *bent* thin
        strip (e.g. a picture-frame annulus going around several corners),
        where the polygon's own oriented-bbox minor extent is dominated by
        the overall bend rather than the material's actual cross-width --
        NOT tied to `spacing` either, since for a strip whose own
        cross-width is smaller than spacing (true for every FinRay ribbon
        here), a spacing-based gap would exceed the strip's own half-width
        and reject every candidate regardless of target count.
    max_rounds: how many doubling batches of Halton candidates to draw
        before giving up.

    Returns an (m, 2) array of accepted points. m may be less than the
    requested target if the polygon is too thin/small to fit that many
    points respecting min_boundary_gap -- this prints a warning rather than
    raising, since a "no fill" ribbon is a valid (if uninteresting) result.
    '''
    exterior = np.asarray(exterior, dtype=float)
    inner_loops = [np.asarray(h, dtype=float) for h in (inner_loops or [])]

    if density_multiplier == 0 and n_target is None:
        return np.empty((0, 2))

    center, u, v, half_u, half_v = oriented_bbox(exterior)
    obb_area = 4.0 * half_u * half_v
    if obb_area < 1e-30:
        return np.empty((0, 2))

    poly_area = _shoelace_area(exterior) - sum(_shoelace_area(h) for h in inner_loops)

    if n_target is None:
        if spacing is None:
            raise ValueError("fill_polygon_with_halton needs `spacing` or `n_target`")
        n_target = max(0, int(round(density_multiplier * poly_area / spacing ** 2)))
    if n_target <= 0:
        return np.empty((0, 2))

    if min_boundary_gap is None:
        perimeter = _perimeter(exterior) + sum(_perimeter(h) for h in inner_loops)
        local_half_width = poly_area / perimeter if perimeter > 1e-30 else half_v
        min_boundary_gap = 0.25 * local_half_width

    path_ext = closed_path(exterior)
    path_inner_loops = [closed_path(h) for h in inner_loops]

    accept_rate_est = max(poly_area / obb_area, 0.05)  # floor avoids an infinite batch
    batch = max(int(np.ceil(n_target / accept_rate_est * 1.5)), 64)

    accepted = []
    start = 1
    have = 0
    for _ in range(max_rounds):
        h01 = halton(batch, dims=2, bases=bases, start=start)
        start += batch

        su = (h01[:, 0] * 2.0 - 1.0) * half_u
        sv = (h01[:, 1] * 2.0 - 1.0) * half_v
        world = center + su[:, None] * u + sv[:, None] * v

        inside = path_ext.contains_points(world)
        for ph in path_inner_loops:
            inside &= ~ph.contains_points(world)
        cand = world[inside]

        if min_boundary_gap > 0 and len(cand):
            cand = cand[_far_from_boundary(cand, exterior, inner_loops, min_boundary_gap)]

        accepted.append(cand)
        have += len(cand)
        if have >= n_target:
            break
        batch *= 2

    pts = np.vstack(accepted) if accepted else np.empty((0, 2))
    if len(pts) > n_target:
        # Halton points are deterministic and low-discrepancy in draw order,
        # so truncating the list (unlike truncating a random point set)
        # still leaves good coverage.
        pts = pts[:n_target]
    elif len(pts) < n_target:
        print(f"[mesh_utils] fill_polygon_with_halton: placed {len(pts)}/{n_target} "
              f"points after {max_rounds} rounds (polygon may be thinner than the "
              f"requested spacing/min_boundary_gap allow); continuing with fewer points.")
    return pts


# --- Boundary-segment recovery (constrained Delaunay) -----------------------

def _orient(a, b, c):
    '''Twice the signed area of triangle abc (>0 if counter-clockwise).'''
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def _in_circumcircle(a, b, c, d):
    '''True if d lies strictly inside the circumcircle of triangle abc.'''
    if _orient(a, b, c) < 0:
        b, c = c, b
    m = np.array([[a[0] - d[0], a[1] - d[1], (a[0] - d[0]) ** 2 + (a[1] - d[1]) ** 2],
                  [b[0] - d[0], b[1] - d[1], (b[0] - d[0]) ** 2 + (b[1] - d[1]) ** 2],
                  [c[0] - d[0], c[1] - d[1], (c[0] - d[0]) ** 2 + (c[1] - d[1]) ** 2]])
    return np.linalg.det(m) > 0


def _triangulate_pseudo_polygon(pts, chain, a, b):
    '''
    Delaunay-triangulates the pseudo-polygon a -> chain... -> b closed by
    segment b-a (Anglada 1997): pick the chain vertex c whose circumcircle
    with a, b holds no other chain vertex, emit (a, b, c), recurse on both
    sides of c.
    '''
    if not chain:
        return []
    c_pos = 0
    for k in range(1, len(chain)):
        if _in_circumcircle(pts[a], pts[b], pts[chain[c_pos]], pts[chain[k]]):
            c_pos = k
    c = chain[c_pos]
    return ([(a, b, c)]
            + _triangulate_pseudo_polygon(pts, chain[:c_pos], a, c)
            + _triangulate_pseudo_polygon(pts, chain[c_pos + 1:], c, b))


def _recover_segment(pts, triangles, p, q):
    '''
    Forces segment p-q into a triangulation that lacks it: removes every
    triangle the open segment crosses, then re-triangulates the cavity on
    each side of p-q. Returns the new triangle list, or the input unchanged
    (plus False) if some vertex lies exactly on p-q (would need splitting).
    '''
    P, Q = pts[p], pts[q]
    crossed, cuts = [], []  # crossed triangles; (t along pq, edge) per crossed edge
    for tri in triangles:
        hit = False
        for u, v in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            if p in (u, v) or q in (u, v):
                continue
            U, V = pts[u], pts[v]
            o1, o2 = _orient(P, Q, U), _orient(P, Q, V)
            o3, o4 = _orient(U, V, P), _orient(U, V, Q)
            if o1 * o2 < 0 and o3 * o4 < 0:
                hit = True
                cuts.append((o3 / (o3 - o4), u, v))
            elif o1 == 0 or o2 == 0:
                w = U if o1 == 0 else V
                s = np.dot(w - P, Q - P) / np.dot(Q - P, Q - P)
                if 0 < s < 1:
                    return triangles, False
        if hit:
            crossed.append(tri)
    # Walk the crossed edges in order from p to q; each has one vertex on
    # either side of p-q, which builds the two cavity chains in order.
    left, right = [], []
    for _t, u, v in sorted(cuts):
        for w in (u, v):
            side = left if _orient(P, Q, pts[w]) > 0 else right
            if w not in side:
                side.append(w)
    crossed_set = {tuple(t) for t in crossed}
    kept = [t for t in triangles if tuple(t) not in crossed_set]
    new = (_triangulate_pseudo_polygon(pts, left, p, q)
           + _triangulate_pseudo_polygon(pts, right, q, p))
    return kept + [tuple(t) for t in new], True


def _recover_loop_segments(pts, triangles, loop_sizes):
    '''
    Ensures every polygon-loop segment is an edge of the triangulation.
    Points 0..sum(loop_sizes)-1 are the loops' vertices, loop after loop,
    each loop closed. Returns (triangles, n_unrecoverable).
    '''
    triangles = [tuple(int(v) for v in t) for t in triangles]
    n_bad = 0
    start = 0
    for n in loop_sizes:
        for k in range(n):
            p, q = start + k, start + (k + 1) % n
            if p == q:
                continue
            present = any(p in t and q in t for t in triangles)
            if not present:
                triangles, ok = _recover_segment(pts, triangles, p, q)
                n_bad += not ok
        start += n
    return triangles, n_bad


# --- Delaunay triangulation -> edge list ------------------------------------

def triangulate_with_holes(exterior, inner_loops, boundary_pts, interior_pts,
                            min_area=1e-14):
    '''
    Triangulates boundary_pts ++ interior_pts with scipy.spatial.Delaunay,
    forces every `exterior`/`inner_loops` segment back in where plain
    Delaunay dropped it (_recover_loop_segments -- on a thin curved ribbon,
    interior points near a long boundary segment otherwise make Delaunay
    connect across that segment, so triangles straddle the wall's edge and
    the centroid test below throws part of the wall area away), then
    discards:
      - any triangle whose centroid falls outside `exterior` or inside any
        `inner_loops` polygon (with every loop segment present, this cleanly
        separates inside from outside triangles)
      - any triangle with area < min_area (numerically-degenerate slivers)

    boundary_pts must be exactly the points of `exterior` followed by the
    points of each `inner_loops[k]`, in that same order and with matching local
    indices 0..nb-1 -- callers use that correspondence to map local
    Delaunay indices back to their own global point indices.
    interior_pts occupy local indices nb..nb+ni-1.

    Returns (all_pts, edges, triangles):
      all_pts: (nb+ni, 2) = vstack(boundary_pts, interior_pts)
      edges: sorted list of unique (i, j) local-index pairs, i < j, one per
          edge of every surviving triangle
      triangles: (nt, 3) array of surviving simplex vertex indices (local)
    '''
    if _Delaunay is None:
        raise ImportError(
            "scipy is required for triangulate_with_holes (scipy.spatial.Delaunay); "
            "see pyIB2d/requirements.txt"
        ) from _SCIPY_IMPORT_ERROR

    boundary_pts = np.asarray(boundary_pts, dtype=float)
    interior_pts = np.asarray(interior_pts, dtype=float) if len(interior_pts) else np.empty((0, 2))
    all_pts = np.vstack([boundary_pts, interior_pts])

    if len(all_pts) < 3:
        return all_pts, [], np.empty((0, 3), dtype=int)

    tri = _Delaunay(all_pts)
    loop_sizes = [len(exterior)] + [len(h) for h in (inner_loops or [])]
    simplices, n_bad = _recover_loop_segments(all_pts, tri.simplices, loop_sizes)
    simplices = np.array(simplices, dtype=int).reshape(-1, 3)
    if n_bad:
        print(f"WARNING (triangulate_with_holes): {n_bad} boundary segment(s) "
              "could not be recovered (a point lies exactly on them)")

    path_ext = closed_path(exterior)
    path_inner_loops = [closed_path(h) for h in (inner_loops or [])]

    edges = set()
    kept = []
    for simplex in simplices:
        p0, p1, p2 = all_pts[simplex]
        centroid = (p0 + p1 + p2) / 3.0
        if not path_ext.contains_point(centroid):
            continue
        if any(ph.contains_point(centroid) for ph in path_inner_loops):
            continue
        area2 = abs((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]))
        if area2 < 2.0 * min_area:
            continue
        kept.append(simplex)
        i, j, k = simplex
        for a, b in ((i, j), (j, k), (k, i)):
            edges.add((min(a, b), max(a, b)))

    triangles = np.array(kept, dtype=int) if kept else np.empty((0, 3), dtype=int)
    return all_pts, sorted(edges), triangles
