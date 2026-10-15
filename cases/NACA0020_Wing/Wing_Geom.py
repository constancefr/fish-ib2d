"""
Builds the Lagrangian structure of the NACA0020 wing from the WaterLily_FSI
example (gvs_wings/Wing_3order.jl) and writes the IB2d input files
(.vertex, .spring, .target, .geo_connect).

Structure model (same method as cases/FinRay/FinRay_Geom.py; see NOTES.md):
the whole airfoil is ONE polygon, meshed once: points every ~h along the
outline, Halton fill inside, one constrained Delaunay triangulation
(build_wing_mesh). Every triangle edge is a spring with the tributary-area
stiffness (common/mesh_utils.tributary_springs), so bending stiffness comes
from the material itself (one side stretches while the other compresses).

  * Reference model: link 2 (0 - 0.06 m from the leading edge) is rigid and,
    because heave/pitch are frozen in the coupled run, fixed in space; link 3
    (0.06 - 0.20 m) is a flexible beam of thickness w(x), E = 1e4 Pa.
  * IB2d: the clamp line x' = 0.06 m is a constrained edge chain of the mesh
    (no triangle crosses it), and every point with x' <= 0.06 m is a target
    point, so the rigid block is held fixed and the flexible part is clamped
    exactly at 0.06 m.
  * The spring network is ~3x softer than E with stiffness_scale = 1 (see
    mesh_utils.tributary_springs); stiffness_scale is measured by
    verify_structure.py.
  * The wing starts pre-bent (init_q, as in the reference): the flexible
    points are mapped onto the bent centerline with plane sections staying
    plane, while every spring keeps the rest length of the straight wing.
  * Everything is per unit depth (span = 1): the coupled WaterLily model scales
    structure (h = 0.1 m) and fluid load (w = 0.1) by the same span, so the
    response is span-independent.

IB2d unit convention (important): the driver spreads `f * ds` with the
*global* ds = Lx / (2 Nx), so a physical nodal force F [N per unit depth]
needs f = F / ds. Every stiffness is therefore computed in physical units and
divided by ds via to_ib() (common/ib_units.py).

The previous model (bending centreline chain + tied skin) is in
Wing_Geom_centerline.py, kept for reference only.
"""

import json
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_COMMON_DIR = os.path.join(os.path.dirname(os.path.dirname(_HERE)), 'common')
for _p in (_HERE, _COMMON_DIR):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import mesh_utils as mu
import wing_case
from ib_units import ib_ds

CHORD = 0.2   # m; the thickness fit below is only valid for this chord


# --------------------------------------------------------------------------
# Airfoil shape
# --------------------------------------------------------------------------

def thickness(xp):
    '''Full airfoil thickness [m] at distance xp [m] from the leading edge.
    Verbatim from Wing_3order.jl (thickness_from_physical_x); zero at the tip.'''
    five_x = 5.0 * np.asarray(xp, dtype=float)
    w = 0.4 * (0.0322 * five_x**3 - 0.0333 * np.cbrt(five_x) - 0.1269 * five_x**2
               + 0.3797 * np.sqrt(five_x) - 0.2518 * five_x)
    return np.maximum(w, 0.0)


def surface_stations(xa, xb, h):
    '''x' [m] of surface stations from xa to xb (both included), spaced ~h in
    *arc length* along the surface (the leading edge is nearly vertical, so
    uniform x' would leave huge gaps there). Both surfaces share them.'''
    u = np.linspace(0.0, 1.0, 200001)
    xp = xa + (xb - xa) * u**2          # refined near xa (the leading edge for xa = 0)
    hw = 0.5 * thickness(xp)
    S = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(xp), np.diff(hw)))])
    n = max(1, int(round(S[-1] / h)))
    xs = np.interp(np.linspace(0.0, S[-1], n + 1), S, xp)
    xs[0], xs[-1] = xa, xb              # exact ends (the clamp station must be exactly xb)
    return xs


def centerline_shape(init_q, L3, s_nodes):
    '''
    Bent centerline of the flexible link for a curvature field
    kappa(X) = (1/L3) * sum_k q_k P_k(2X-1), X = s/L3 (same Legendre basis and
    1/L scaling as the GVS model), sampled at arc lengths s_nodes.
    Returns (x, y, theta) relative to the clamp, tangent = +x at s = 0.
    '''
    fine = np.linspace(0.0, L3, 20001)
    kappa = np.polynomial.legendre.legval(2.0 * fine / L3 - 1.0, init_q) / L3
    theta = np.concatenate([[0.0], np.cumsum(0.5 * (kappa[1:] + kappa[:-1]) * np.diff(fine))])
    xf = np.concatenate([[0.0], np.cumsum(0.5 * (np.cos(theta[1:]) + np.cos(theta[:-1])) * np.diff(fine))])
    yf = np.concatenate([[0.0], np.cumsum(0.5 * (np.sin(theta[1:]) + np.sin(theta[:-1])) * np.diff(fine))])
    return np.interp(s_nodes, fine, xf), np.interp(s_nodes, fine, yf), np.interp(s_nodes, fine, theta)


def _wing_outline(h, rigid_len):
    '''
    The airfoil in the wing frame (x' from the leading edge along the chord,
    y' up): one closed polygon LE -> upper surface -> TE -> lower surface
    (first vertex not repeated), points every ~h in arc length, with a vertex
    on each surface exactly at the clamp station x' = rigid_len.

    Returns (outer, clamp_pts, clamp_ends):
        outer: (n, 2) outline polygon
        clamp_pts: (m, 2) points strictly inside the clamp line (upper ->
            lower at x' = rigid_len), spaced ~h
        clamp_ends: (upper, lower) indices of the clamp line's ends in outer
    '''
    xs = np.concatenate([surface_stations(0.0, rigid_len, h)[:-1],
                         surface_stations(rigid_len, CHORD, h)])
    hw = 0.5 * thickness(xs)
    n = len(xs) - 1                                   # xs[0] = LE, xs[n] = TE
    upper = np.column_stack([xs[1:n], hw[1:n]])
    lower = np.column_stack([xs[n - 1:0:-1], -hw[n - 1:0:-1]])
    outer = np.vstack([[0.0, 0.0], upper, [CHORD, 0.0], lower])

    k = int(np.argmin(np.abs(xs - rigid_len)))        # station index of the clamp
    i_up = k                                          # outer index of upper point k (LE is 0)
    i_lo = len(outer) - k                             # outer index of lower point k
    half = hw[k]
    m = max(1, int(round(2.0 * half / h)))
    ys = np.linspace(half, -half, m + 1)[1:-1]
    clamp_pts = np.column_stack([np.full(len(ys), rigid_len), ys])
    return outer, clamp_pts, (i_up, i_lo)


def build_wing_mesh(fill_spacing, fill_density_multiplier, rigid_len):
    '''
    Meshes the whole airfoil at once: points every ~h along the outline and
    the clamp line, Halton fill inside, one constrained triangulation with the
    clamp line as a chain of extra edges. h = fill_spacing / sqrt(multiplier)
    is the fill's mean spacing, so edge and interior points are equally dense.

    Halton points keep min_boundary_gap = h/2 from the outline and the clamp
    line. (mu's default, 0.25 * area / perimeter, would be ~3.3 mm here: more
    than ds and than the trailing edge's thickness, leaving a fill-free band
    along the skin.)

    Returns (pts, triangles, info): pts (N, 2) in the wing frame (outline
    first, then the clamp line, then the fill); info holds n_outline,
    clamp_idx (the clamp line's indices, upper end to lower end), te_index,
    outer (the outline polygon) and h.
    '''
    h = fill_spacing / np.sqrt(max(fill_density_multiplier, 1.0))
    outer, clamp_pts, (i_up, i_lo) = _wing_outline(h, rigid_len)
    fill = mu.fill_polygon_with_halton(outer, fill_spacing,
                                       density_multiplier=fill_density_multiplier,
                                       min_boundary_gap=0.5 * h,
                                       keep_out_segments=[(outer[i_up], outer[i_lo])])
    nb = len(outer)
    clamp_idx = [i_up] + list(range(nb, nb + len(clamp_pts))) + [i_lo]
    pts, _, triangles = mu.triangulate_with_holes(
        outer, [], outer, np.vstack([clamp_pts, fill]),
        extra_segments=list(zip(clamp_idx[:-1], clamp_idx[1:])))
    te_index = int(np.argmin(np.hypot(pts[:, 0] - CHORD, pts[:, 1])))
    return pts, triangles, dict(n_outline=nb, clamp_idx=clamp_idx, te_index=te_index,
                                outer=outer, h=h)


def bend(pts, init_q, rigid_len):
    '''Initial (pre-bent) coordinates: every point of the flexible link at
    (x', y') goes to c(s) + y' n(s), s = x' - rigid_len, with c, n the bent
    centerline and its normal (plane sections stay plane and normal to the
    centerline, as in the reference's beam kinematics).'''
    out = pts.copy()
    flex = pts[:, 0] > rigid_len
    if not flex.any() or not any(abs(q) > 0 for q in init_q):
        return out
    cx, cy, th = centerline_shape(init_q, CHORD - rigid_len, pts[flex, 0] - rigid_len)
    yp = pts[flex, 1]
    out[flex, 0] = rigid_len + cx - yp * np.sin(th)
    out[flex, 1] = cy + yp * np.cos(th)
    return out


# --------------------------------------------------------------------------
# Structure builder (importable by verify_structure.py)
# --------------------------------------------------------------------------

def build_wing(params=None, Nx=None, Lx=None, **overrides):
    '''
    Returns a dict with the vertex coordinates (x, y: bent initial state;
    xr, yr: straight rest state, both in world coordinates), the springs
    (connections, edge_k, in IB2d units; all_connections = every mesh edge,
    including those inside the clamped block), the target points and their
    k_Target, and mesh bookkeeping. `overrides` replace geom keys of
    case_params.json (e.g. stiffness_scale=1.0 for calibration).
    '''
    g = dict((params or wing_case.load_case_params()).get('geom', {}))
    g.update(overrides)
    Nx = Nx if Nx is not None else int(wing_case.read_input2d('Nx'))
    Lx = Lx if Lx is not None else wing_case.read_input2d('Lx')
    ds = ib_ds(Lx, Nx)            # IB2d's global Lagrangian spacing (hard-coded in the driver)

    E = g.get('E_flex', 1.0e4)
    rigid_len = g.get('rigid_length', 0.06)
    init_q = list(g.get('init_q', [0, 0, 0, 0]))
    mult = g.get('fill_density_multiplier', 1.0)
    fill_spacing = g.get('fill_spacing') or ds        # the json holds null -> ds
    stiffness_scale = g.get('stiffness_scale', 1.0)
    target_factor = g.get('target_stiffness_factor', 15.0)
    px, py = g.get('pitch_axis_world', [0.4, 0.6])
    x_LE = px - g.get('pitch_axis_from_LE', 0.03)
    y0 = py

    t_start = time.perf_counter()
    pts, triangles, info = build_wing_mesh(fill_spacing, mult, rigid_len)
    xr, yr = pts[:, 0], pts[:, 1]
    connections, edge_k = mu.tributary_springs(xr, yr, triangles, E, ds, stiffness_scale)

    # clamp: every point of the rigid block, clamp line included
    target_idx = np.where(xr <= rigid_len + 1e-12)[0]
    # a spring between two clamped points never stretches: drop it (IB2d's
    # spring force is a Python loop, ~6 us per spring per call; this removes
    # about a third of them)
    rigid = set(target_idx.tolist())
    all_connections = connections        # still drawn: .geo_connect and the preview
    connections = [c for c in connections if not (c[0] in rigid and c[1] in rigid)]
    edge_k = {(i, j): edge_k[(i, j)] for i, j, _ in connections}
    k_Target = target_factor * mu.max_incident_stiffness(target_idx.tolist(), edge_k)
    bent = bend(pts, init_q, rigid_len)
    build_time = time.perf_counter() - t_start

    return dict(x=bent[:, 0] + x_LE, y=bent[:, 1] + y0, xr=xr + x_LE, yr=yr + y0,
                triangles=triangles, connections=connections, edge_k=edge_k,
                all_connections=all_connections,
                target_idx=target_idx, k_Target=k_Target, ds=ds, fill_spacing=fill_spacing,
                mult=mult, E=E, stiffness_scale=stiffness_scale, rigid_len=rigid_len,
                x_LE=x_LE, y0=y0, build_time=build_time, **info)


def mesh_info(w):
    '''Mesh metrics for the run record (run_case.py copies them into run_info.json).'''
    area = mu._shoelace_area(w['outer'])
    x, y = w['xr'][w['triangles']], w['yr'][w['triangles']]
    mesh_area = 0.5 * np.abs((x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0])
                             - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0])).sum()
    k = np.array(list(w['edge_k'].values()))
    n_edge = w['n_outline'] + len(w['clamp_idx']) - 2
    return dict(
        ds=w['ds'], fill_spacing=w['fill_spacing'], fill_density_multiplier=w['mult'],
        edge_spacing=w['h'], n_edge_points=n_edge,
        n_fill_points_requested=int(round(w['mult'] * area / w['fill_spacing'] ** 2)),
        n_fill_points=len(w['xr']) - n_edge,
        n_points=len(w['xr']), n_triangles=len(w['triangles']), n_springs=len(w['connections']),
        stiffness_scale=w['stiffness_scale'], k_median=float(np.median(k)), k_max=float(k.max()),
        n_target_points=len(w['target_idx']), n_clamp_points=len(w['clamp_idx']),
        k_Target=w['k_Target'], te_index=w['te_index'],
        wing_area_m2=area, mesh_area_fraction=mesh_area / area,
        build_time_s=round(w['build_time'], 3),
    )


# --------------------------------------------------------------------------
# File writers (same formats as FinRay_Geom.py)
# --------------------------------------------------------------------------

def write_files(w, struct_name='wing'):
    with open(struct_name + '.vertex', 'w') as f:
        f.write('%d\n' % len(w['x']))
        for xi, yi in zip(w['x'], w['y']):
            f.write('%1.16e %1.16e\n' % (xi, yi))

    # rest length = the STRAIGHT wing's edge length, so a pre-bent start is
    # the initial elastic state
    with open(struct_name + '.spring', 'w') as f:
        f.write('%d\n' % len(w['connections']))
        for i, j, _ in w['connections']:
            rl = float(np.hypot(w['xr'][j] - w['xr'][i], w['yr'][j] - w['yr'][i]))
            f.write('%d %d %1.16e %1.16e %1.16e\n' % (i, j, w['edge_k'][(i, j)], rl, 1.0))

    with open(struct_name + '.target', 'w') as f:
        f.write('%d\n' % len(w['target_idx']))
        for i in w['target_idx']:
            f.write('%d %1.16e\n' % (i, w['k_Target']))

    with open(struct_name + '.geo_connect', 'w') as f:   # every mesh edge (drawing only)
        for i, j, _ in w['all_connections']:
            f.write('%d %d\n%d %d\n' % (i, j, j, i))


def plot_preview(w, path='geometry_preview.png'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 1, figsize=(12, 9))
    nf = w['n_outline'] + len(w['clamp_idx']) - 2      # first fill point
    for a, xx, yy, title in ((ax[0], w['xr'], w['yr'], 'rest (straight)'),
                             (ax[1], w['x'], w['y'], 'initial (pre-bent)')):
        for i, j, kind in w['all_connections']:
            a.plot([xx[i], xx[j]], [yy[i], yy[j]], color='b' if kind == 'boundary' else (0.6, 0.6, 0.6, 0.6),
                   lw=0.75 if kind == 'boundary' else 0.4)
        a.plot(xx[nf:], yy[nf:], 'k.', ms=2)
        a.plot(xx[w['target_idx']], yy[w['target_idx']], 'r.', ms=4, label='target (clamped)')
        a.plot(xx[w['clamp_idx']], yy[w['clamp_idx']], 'm-', lw=1.5, label='clamp line')
        a.plot(xx[w['te_index']], yy[w['te_index']], 'g*', ms=10, label='trailing edge')
        a.set_aspect('equal'); a.set_title('%s: %d points, %d springs, x%g'
                                           % (title, len(xx), len(w['connections']), w['mult']))
    ax[0].legend(loc='upper right', fontsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=150)


def main():
    w = build_wing()
    write_files(w)
    info = mesh_info(w)
    with open('mesh_info.json', 'w') as f:
        json.dump(info, f, indent=2)
    print('%d points (%d edge, %d fill), %d springs, %d target points, k_Target=%.3e, '
          'mesh area %.4f%% (ds=%.5f, h=%.5f)'
          % (info['n_points'], info['n_edge_points'], info['n_fill_points'], info['n_springs'],
             info['n_target_points'], info['k_Target'], 100 * info['mesh_area_fraction'],
             info['ds'], info['edge_spacing']))
    plot_preview(w)


if __name__ == '__main__':
    main()
