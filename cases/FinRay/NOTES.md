In FinRay_Geom.py -----------

================================================================================
SUMMARY (high-level overview of this round of work -- see the matching
detailed section below, by heading name, for exact numbers and how each
finding was verified)
================================================================================

Geometry & meshing:
  - Added Halton-sampled interior "fill" points + Delaunay triangulation to
    each ribbon (outline + 4 rays): _mesh_utils.py (Halton generator,
    oriented-bbox sampler, polygon fill, triangulation) plus
    case_params.json["geom"]'s ribbon_fill_spacing/
    ribbon_fill_density_multiplier/ribbon_fill_stiffness_scale. See
    "ribbon_fill_density_multiplier -> actual point count" below for the
    formula and concrete numbers (multiplier=1 is only 6 points total at
    this example's defaults -- the ribbons are thin relative to ds).
  - Replaced each ray's old "spring to the nearest existing outline point"
    'attach' connection with a true weld (shared vertex) -- see
    "build_Tail_Beams" below. Removed the 'attach' kind and its
    k_Attach=50*k_Spring. This also incidentally fixed a NaN-instability
    (k_Attach was the model's stiffest element; see
    "build_Tail_Ribbon_Connections" below).
  - Found + fixed a rail-desync bug the weld introduced: 'rung'/'diag'
    springs connecting far-apart points instead of ones actually across
    the ribbon (visible as spurious diagonals cutting through the tail's
    open "windows"). See "build_Tail_Beams" below.
  - Fixed boundary point-count rounding (_n_points_for_length): confused
    "number of points" with "number of segments" for short open chains
    (base, rays), leaving some segments up to 214% of ds. Now 70-124%
    across the whole structure.
  - Found, NOT fully solved: 'fill' triangles can invert under large
    actuation-driven bending -- small, localized, self-correcting, and a
    fundamental limitation of pure spring lattices (nothing prevents
    inversion, only resists it) rather than a bug. Worst at the base. See
    "'fill' triangle inversion under large deformation" below.

Stiffness model:
  - (2026-10-05, SUPERSEDES the rail/rung/diag formulas below) Every
    spring now uses ONE formula, k = scale*E*depth*A_tributary/L^2, from
    one triangulation per ribbon; separate 'rung'/'diag' springs are gone.
    The old scheme counted the wall's material ~1.75x (density x1) to
    ~2.3x (x64), growing with density; now exactly 1x at every density.
    ribbon_fill_stiffness_scale now scales ALL springs. Also fixed two
    triangulation bugs that silently dropped wall area. See "single
    stiffness formula" below.
  - Spring stiffness is now computed per-spring from its own rest length
    (k = E*A/L_own, standard truss-element formula), not one constant per
    kind -- needed once weld-spliced springs could be much shorter than
    the nominal ds. See "print_Lagrangian_Springs" below.
  - ribbon_fill_stiffness_scale=1.0 (the formula's literal value) caused
    NaN before the weld fix; clean after it, at the densities tested so
    far. The fill-spring formula also has a hard ceiling independent of
    scale: a central-force triangular lattice can't reproduce silicone's
    near-incompressible Poisson's ratio (pinned at ~1/3). A fluid-free
    Timoshenko-beam calibration was recommended but NOT YET BUILT -- see
    "ribbon_fill_stiffness_scale -- what it represents" below.

Actuation:
  - Fixed a case-sensitivity bug in update_Target_Point_Positions.py
    (read lowercase 'frequency' etc. while case_params.json uses uppercase
    FREQUENCY) that silently no-op'd every actuation override.
  - Each base corner's actuation is now a small rigid PATCH of nearby
    boundary points (both rails + anything within actuation_patch_radius),
    not a single point. Radius defaults to wall_thickness, NOT ds -- a
    ds-sized radius was found to reach into the adjacent edge and
    destabilize the sim. See "_actuation_patch_indices" below.
  - (2026-10-06) k_Target is now 15 x the largest total spring stiffness
    attached to any patch point (geom key target_stiffness_factor),
    replacing 1000*k_Spring in FinRay_Geom.py x TARGET_STIFFNESS_SCALE=0.2
    at runtime (now removed). .beam file, k_Spring, k_Beam and the beam
    helpers removed. See "beams removed / k_Target" below.

Checked and found NOT to be a problem:
  - input2d's beams=0 predates all of this work (set since the example's
    very first commit) -- not a consequence of the current mesh. Verified
    directly (rest-vs-late turning-angle comparison) that the rung/diag
    truss alone already prevents the zig-zag artifact beams exist to
    guard against, with or without the fill mesh. See "beams flag
    (input2d) / zig-zag check" below.

Infrastructure (outside FinRay_Geom.py):
  - run_case.py: BASE_CASE now loads from case_params.json (not a
    separate hardcoded dict); 2D sweep over density_multiplier x
    stiffness_scale; automatic NaN-divergence detection per case.
  - make_movie.py: new `--all` mode sweeps every results/ subfolder in
    one VisIt session, with per-case database cleanup and failure
    isolation (one bad case doesn't abort the rest).

Bugs found from incomplete prior edits (FOUND + FIXED):
  - An in-progress holes->inner_loops rename was completed across
    FinRay_Geom.py but not _mesh_utils.py, so fill_polygon_with_halton()
    raised TypeError on every single run. Fixed.
  - _Build_Tail_Geometry_Ribbon's `fill_spacing is None -> ds` fallback
    had been deleted even though its own docstring still promised it,
    breaking any caller (e.g. give_Me_Immersed_Boundary_Geometry) that
    passes fill_spacing=None directly. Restored.
  See "Incomplete-rename / dropped-fallback bugs" below for both.

================================================================================
DETAILED NOTES
================================================================================

beams removed / k_Target (2026-10-06):
    WHAT CHANGED. FinRay_Geom.py no longer writes finray.beam, and
    k_Spring, k_Beam, build_Tail_Ribbon_Beams, print_Lagrangian_Beams,
    spring_Stiffness_From_Material and beam_Stiffness_From_Material are
    gone (all in git history). k_Target is now
        k_Target = target_stiffness_factor * max over patch points of
                   (sum of k of every spring attached to that point)
    with target_stiffness_factor = 15 (case_params.json "geom").
    TARGET_STIFFNESS_SCALE is removed from update_Target_Point_Positions.py
    / case_params.json "actuation"; the .target value is used as-is (still
    ramped in over STIFFNESS_RAMP_TIME).

    WHY. input2d has beams=0, so .beam was never read, and the zig-zag
    check ("beams flag" below) showed it isn't needed. Once beams went,
    k_Spring (a formula no spring uses since "single stiffness formula")
    only survived to size k_Target. A fixed absolute k_Target was rejected
    too: what matters is k_Target RELATIVE to the springs it must overpower
    (too low -> patch lags the prescribed motion; too high -> explicit
    time-stepping NaN), and those move with E_material, stiffness_scale,
    mesh density and ds. Global max(edge_k) was also rejected as the
    reference -- it's set by whichever sliver triangle the Halton fill
    happens to make: 4.9e4 (x0, x1), 1.0e5 (x4), 1.5e5 (x8), 9.4e4 (x16),
    6.9e4 (x32), 6.5e4 (x64). The incident-stiffness reference is steady:
    3.5e4 (x0, x1), 5.8e4-6.9e4 (x4 .. x64).

    VALUE. Old effective k_Target = 0.2 * 1000 * k_Spring = 9.38e5 (E=0.74e6
    defaults), i.e. 13.5-16x the incident stiffness at x4..x64 and 27x at
    x0/x1. Factor 15 reproduces the x4..x64 behaviour; at x0/x1 k_Target
    roughly halves (5.2e5) -- re-check patch tracking there if it matters.

single stiffness formula (build_Tail_Ribbon_Connections, 2026-10-05):
    WHAT CHANGED. Every spring is an edge of its ribbon's triangulation
    (rail points + interior Halton points; triangulated even at density
    multiplier 0) and gets
        k = stiffness_scale * E * extrude_depth * A_tributary / L**2,
    A_tributary = 1/3 of each triangle sharing the edge, summed over every
    ribbon the edge belongs to. Kinds 'frame'/'ray' (rails) and 'fill'
    (everything else) are kept for plotting only. 'rung'/'diag' and
    shear_Stiffness_From_Material are removed; EA_by_kind is gone from
    print_Lagrangian_Springs. ribbon_fill_stiffness_scale keeps its name
    but now multiplies every spring (the single calibration knob alpha).
    k_Spring = E*(t/2)*depth/ds survived only to size k_Beam/k_Target
    (since removed -- see "beams removed / k_Target").

    WHY. "Material budget" check: sum over springs of k*L^2/(E*depth)
    equals the wall area each spring set represents (exact for the
    tributary formula). Old scheme, as % of the true wall area:
    rails ~100% (each rail = half the wall) + rung/diag ~50% + fill 23%
    (x1) .. 82% (x64) = 175% .. 234% -- material double-counted, and the
    fin got stiffer as density rose, which would confound any density
    study. New scheme: 100% at x0, x1, x4, x16, x64, x256 (rails 36% ->
    5% as their edge strip thins, fill the rest). Rails now get roughly
    half the k of a similar-length interior edge (triangles on one side
    only). 100% is necessary, not sufficient: actual bending stiffness
    still has to be checked against beam theory / FE (calibration step).

    TRIANGULATION BUGS FOUND + FIXED (in _mesh_utils.py) while checking:
    1. matplotlib Path(loop, closed=True) IGNORES the last vertex (treats it
       as the CLOSEPOLY placeholder). Every inside/outside test dropped one
       corner per loop: triangles near it were discarded (one ray ribbon
       lost ~70% of its area at x16) and Halton points near it rejected.
       New helper mu.closed_path() appends the first vertex; used in
       fill_polygon_with_halton, triangulate_with_holes, and the rail
       orientation test in _Build_Tail_Geometry_Ribbon.
    2. Plain (unconstrained) Delaunay on the thin, curved frame ribbon
       drops boundary segments: a 15 mm rail segment with interior points
       within its diametral circle gets crossed by edges instead, those
       triangles straddle the wall edge, and the centroid filter discards
       them (13-16 of 70 frame segments missing, ~9% of frame area lost,
       at x16-x64). triangulate_with_holes now recovers every loop segment
       (_recover_loop_segments: remove crossed triangles, re-triangulate
       each side's cavity, Anglada 1997) -- i.e. constrained Delaunay
       without a new dependency (`triangle` has no py3.14 wheel). After
       both fixes: 0 missing segments and 100% area coverage, x0..x256.

    KNOWN, NOT FIXED: rail points stay ~ds (~15 mm) apart at every density
    while interior spacing shrinks (~2 mm at x64), so edge triangles get
    long and thin at high density. Spring k range ~150 .. ~1e5 (similar
    max to before). Consider refining the rails with the fill density
    for the density study.

_actuation_patch_indices / actuation_patch_radius (base actuation as a
patch, not a single point):
    Each base corner used to have its actuation displacement applied to a
    SINGLE Lagrangian point (and, because of the update_Target_Point_
    Positions.py argmax/argmin row lookup, actually only ONE of that
    corner's 2 rail points -- the other sat completely fixed, never
    updated). Concentrating the whole prescribed motion onto one point
    while its immediate ~wall_thickness-away neighbor stayed rigid was
    producing a visible local twist right at the corner, distinct from the
    real Fin Ray shear/twist that alternately driving the two corners is
    supposed to produce (see update_Target_Point_Positions.py's own
    docstring -- the alternating-corner scheme itself is intentional and
    unchanged; this only addresses how each corner's own patch is driven).

    Fixed by _actuation_patch_indices: expands each corner into every
    Lagrangian point within actuation_patch_radius of that corner's true
    position, restricted to BOUNDARY (rail) points only (ribbonPolys'
    interior_idx points are excluded). update_Target_Point_Positions.py
    then classifies every target row into the "top" or "bottom" patch by
    baseline y vs. the midline, and moves every point in a patch by the
    identical rigid displacement -- approximating a real actuator's finite
    contact area. Verified directly: at a given phase, every row in a
    patch gets bit-identical dx (np.allclose), confirming rigid (not
    independent) motion.

    STABILITY WARNING (confirmed, not hypothetical): the natural first
    choice, actuation_patch_radius defaulting to ~ds, is WRONG for this
    geometry and was rejected. At ds-scale radius, the patch pulled in the
    nearest point on the ADJACENT top/bottom slant edge (~13-14mm away)
    rather than more of the base edge itself (whose own next point is
    ~17mm away, i.e. *further* than the adjacent edge's point) -- ds
    happens to be roughly comparable to both distances for this
    geometry's proportions, so nothing guarantees a ds-sized radius stays
    within the same edge. Forcing that adjacent-edge point into rigid
    lockstep with the corner (a location it doesn't physically belong to)
    measurably destabilized the simulation: reproduced going fully NaN
    (case: E_material=0.74e6, dt=1e-3, ribbon_fill_density_multiplier=32,
    ribbon_fill_stiffness_scale=1.0, FREQUENCY=0.1, Tfinal=5.0) with an
    8-point patch (ds-radius, boundary-only), while the exact same case
    stayed clean through the full run with a 4-point patch (both rails of
    the true corner only, no extra points).

    Fixed by defaulting actuation_patch_radius to 1.0*wall_thickness
    instead of ds -- ties the radius to the corner's own actual physical
    scale (the 2 rails are only ~wall_thickness apart) rather than the
    mesh resolution, which for this geometry keeps it safely under the
    ~13mm distance to the nearest unrelated point on an adjacent edge.
    Re-verified stable (same Tfinal=5.0 case) at this default.

    Real, NOT-yet-solved limitation: this is a Euclidean-distance patch,
    not an edge-aware (arc-length-along-the-same-rail-chain) one, so
    future geometry changes (different wall_thickness/rayPositions/L/W) could in
    principle put another edge's point closer than intended again --
    there's no general guarantee a given actuation_patch_radius stays
    within the base edge for arbitrary geometry. If you change
    actuation_patch_radius (or the base geometry) from these defaults,
    re-check stability the same way this was found: run to Tfinal and
    check run_case.py's own divergence warning (or grep the last
    viz_IB2d/lagPtsConnect.*.vtk frame for "nan").

_n_points_for_length (boundary density non-uniformity, FOUND + FIXED):
    The base, top/bottom slant edges, and each of the 4 rays used to each
    independently compute their own point count as
    n_points = max(round(length/ds), 2). This looks reasonable but is
    wrong for an OPEN chain: n points give n-1 SEGMENTS, and rounding the
    POINT count directly silently drops that -1. Negligible when there are
    many segments (the slant edges, 8-9 points, ended up ~110% of ds --
    fine), but severe when there are few (the base and rays, 2-3 points):
    measured actual segment lengths (at this example's default ds) ranged
    up to 214% of ds (ray at d=0.047 -- 2 points forced by rounding down
    from 2.14, giving exactly 1 segment spanning its FULL length), with
    the base itself at 163% of ds. Boundary density was NOT actually
    uniform across the structure, worst specifically at the base and the
    inner rays -- exactly the regions carrying the most actuation-driven
    deformation (see the triangle-inversion note below).

    Fixed by rounding to the nearest number of SEGMENTS first
    (n_segments = max(round(length/ds), 1)), then deriving point count
    from that (n_points = n_segments + 1) -- applied uniformly to the
    base, both slant edges, and every ray via the new
    _n_points_for_length() helper. Result: base segments now 94-124% of
    ds (was 163%), rays now 70-96% of ds (was up to 214%), slant edges
    ~98% of ds (was ~110%, also improved). Verified: no degenerate
    springs at any tested ribbon_fill_density_multiplier, and re-ran the
    triangle-inversion check below (same case: mult=32, scale=1.0,
    FREQUENCY=0.1, Tfinal=5.0) -- inversions are delayed and reduced but
    NOT eliminated (outline: first inversion t=0.44s->0.62s, peak
    5->4/378 triangles; ray0: t=0.9s->1.58s, peak 4->1/41 triangles).
    Consistent with the inversion note's own conclusion: this fix removes
    unnecessary sparseness (one contributing factor), but doesn't address
    the deeper limitation that a pure spring network has no mechanism to
    actually PREVENT inversion, only resist it proportionally to local
    stiffness.

'fill' triangle inversion under large deformation (KNOWN, not a bug):
    During simulation (not the static geometry_preview.png), some 'fill'
    triangles in the outline and ray meshes can invert (flip orientation)
    under large actuation-driven bending -- visible as an interior Halton
    point ending up geometrically on the wrong side of a nearby boundary
    edge. Confirmed by tracking each fill triangle's signed area over a
    full simulation run (case: E_material=0.74e6, dt=1e-3,
    ribbon_fill_density_multiplier=32, ribbon_fill_stiffness_scale=1.0,
    FREQUENCY=0.1, Tfinal=5.0) and comparing sign to its rest
    configuration. Before the _n_points_for_length fix above: outline's
    first inversion at t=0.44s (peak 5/372 triangles concurrently
    inverted), ray0's (the widest ray) at t=0.9s (peak 4/38). Both are
    small, localized, and self-correcting -- 'fill' springs are linear, so
    an inverted triangle isn't permanently stuck, it gets pulled back.

    Root cause: a pure mass-spring lattice has no explicit constraint
    against a triangle folding over -- springs only RESIST it in
    proportion to local stiffness and mesh density, they don't PREVENT it
    (unlike e.g. a co-rotational FEM formulation or an explicit
    area-preservation/inversion-barrier term, neither of which this model
    has). The base is hit first/worst because it's simultaneously (a) the
    most forced part of the structure -- actuation drives it directly via
    the target points -- and (b) was also the most sparsely-resolved
    boundary in the whole structure (see _n_points_for_length above),
    meaning larger, less-braced triangles there.

    _n_points_for_length measurably helps (delays onset, reduces peak
    concurrent inversions -- see that note for numbers) by removing the
    unnecessary sparseness, but does not eliminate inversion entirely,
    since the underlying "nothing prevents it" limitation is still there.
    If it becomes more than cosmetic, further options (not yet
    implemented): push ribbon_fill_density_multiplier/
    ribbon_fill_stiffness_scale higher specifically in the base region, or
    a more invasive fix like an explicit anti-inversion penalty term.

print_Lagrangian_Springs (length-aware stiffness):
    Each 'outline'/'ray'/'rung'/'diag' spring's stiffness is computed
    per-spring as k = EA_by_kind[kind] / ds_Rest (ds_Rest = that spring's
    OWN measured rest length), not a single k_Spring/k_Rung constant
    shared by every spring of that kind. This is the standard truss-element
    formula: for a rod of Young's modulus E, cross-section A, rest length
    L, under strain eps=dL/L, force = E*A*eps = (E*A/L)*dL, so
    k = E*A/L -- a SHORTER piece of the same material is stiffer (less
    material to stretch across), a LONGER piece is softer.

    Why this matters here specifically: FinRay_Geom() used to compute
    k_Spring = E*A/ds ONCE using the nominal ds, and hand that same number
    to every 'outline'/'ray' spring regardless of its own actual length.
    Most rail segments ARE close to ds, so this was a fine approximation
    for them -- but weld_to_inner_rail's junction splicing creates some
    springs as short as ~6% of ds. Under the old constant-k scheme those
    were silently ~17x too soft for their length (k should scale as 1/L,
    but they got the same k as a full-ds spring). Fixed by having
    FinRay_Geom() pass EA_by_kind (E*A or G*A products, NOT yet divided by
    a length -- recovered from the existing k_Spring/k_Rung via
    E*A = k_Spring*ds, so spring_Stiffness_From_Material/
    shear_Stiffness_From_Material didn't need to change) and computing
    k = EA_by_kind[kind] / ds_Rest inside print_Lagrangian_Springs, per
    spring. k_Spring/k_Rung themselves are still used as-is for k_Beam/
    k_Target, which want one representative magnitude, not a per-spring
    value. Verified: springs near 100% of ds keep ~the same k as before;
    a 6%-of-ds spring's k rose ~17x; a spring at ~190% of ds dropped to
    ~half -- all consistent with k=EA/L. Re-ran the fluid dynamics check
    (Tfinal=5.0, ribbon_fill_density_multiplier=1.0,
    ribbon_fill_stiffness_scale=0.1) after this change: still clean, no
    NaN, despite several springs getting substantially stiffer.

beam_Stiffness_From_Material:
    IB2d's torsional spring energy is defined via a signed-area
    form of "curvature" rather than a normalized curvature. Might have to
    verify that the scaling of k_Beam with E and I is correct for IB2d's
    discrete implementation.

    Recommended calibration before trusting comparisons between designs:
    build an isolated single strip (no ribs, no fluid), fix one end with
    stiff target points, apply a known force at the free tip, and check
    the resulting deflection against the analytical cantilever formula
        delta = F * L^3 / (3 * E * I)
    Adjust k_Beam until simulated and analytical deflection agree.

build_Tail_Beams:
    This does not add a beam spanning the
    ray-to-outline junction itself, so bending moment isn't transmitted
    across that connection (only the axial attachment spring force is
    via build_Tail_Connections). The junction currently behaves like a hinge
    instead of a fused connection.

    RESOLVED for the current ribbon model: _Build_Tail_Geometry_Ribbon's
    weld_to_inner_rail welds each ray rail's two ends directly into the
    outline's own inner-rail point sequence (a shared vertex, found via
    _mesh_utils.project_point_to_polyline -- the true point where the ray's
    edge meets the outline, not just the nearest existing outline point),
    instead of adding a separate 'attach' spring to a nearby-but-different
    point. The 'attach' kind and its k_Attach=50*k_Spring stiffness no
    longer exist. This still isn't a full moment-transmitting weld (no beam
    spans across the junction, so relative rotation between the ray and the
    outline right at that shared point is still only resisted by whatever
    rung/diag/fill springs happen to reach it, not by a beam), but the
    standoff gap + isolated one-spring joint is gone.

    BUG FOUND + FIXED during this change: weld_to_inner_rail originally
    only spliced new junction points into the outline's INNER rail
    (hole_pts/hole_idx), leaving the OUTER rail (ext_pts/ext_idx) at its
    original length. add_rungs_and_diagonals (in
    build_Tail_Ribbon_Connections) pairs the two rails by matching array
    position k, assuming len(railA) == len(railB) and that railA[k] sits
    directly across the ribbon from railB[k] -- true before any splicing
    (both rails come from the same centerline via _offset_polyline), but
    silently violated from the first inner-rail-only splice onward: railB
    now had extra points the array-index pairing didn't account for, so
    'rung'/'diag' springs past that point connected increasingly
    mismatched (physically far-apart) pairs -- visible as long spurious
    diagonal lines cutting straight across the tail's open "windows" in
    geometry_preview.png. Fixed by having weld_to_inner_rail splice a
    matching "partner" point into the outer rail too, at the same segment
    index and parametric position t (accurate given wall_thickness is tiny
    relative to the outline's feature sizes -- same justification used for
    the Delaunay centroid-filter approximation elsewhere in this file) --
    the partner point has no role of its own, it exists purely to keep
    both rails the same length with the same per-position correspondence.

beams flag (input2d) / zig-zag check (INVESTIGATED, left as-is):
    input2d's beams=0 has been set since this example's very first commit
    (checked via `git log`/`git show` on input2d) -- predates the ribbon
    model, Halton fill, and everything else in this file. It is NOT a
    response to "the current mesh tessellation" and was never toggled as
    part of any of this work.

    What the (never-read, and since 2026-10-06 no longer written) .beam file represents:
    per build_Tail_Ribbon_Beams's own design, a weak (k_Beam =
    beam_reg_fraction * k_Spring = 5% of k_Spring, NOT derived from a real
    EI) per-rail anti-zigzag regularization -- explicitly NOT the source of
    the tail's bending stiffness (that comes from the rail truss itself).
    So enabling beams is not an "accuracy vs. real silicone" lever
    regardless of whether it's needed -- it only guards against a specific
    IBM discretization artifact (a point chain buckling into a sawtooth
    under local compression instead of bending smoothly).

    Checked directly whether that artifact is actually happening, rather
    than assuming: tracked the SIGNED TURNING ANGLE at every point along
    the outline-frame rail and the ray0 rail, comparing rest (t=0) to late
    in a 15s run (t=15s; case: E_material=0.74e6, dt=1e-3, FREQUENCY=0.1,
    AMPLITUDE=0.01), in both the sparsest-mesh regime
    (ribbon_fill_density_multiplier=0 -- bare rung/diag truss, zero fill
    bracing, the case most likely to need beams) and the densest tested
    (multiplier=32). Result in every case: turning angles stayed smooth
    and small (max change from rest ~1.4-1.6 degrees, excluding the real
    geometric corners at the base/tip), with NO alternating-sign sawtooth
    pattern -- i.e. no zig-zag, with or without the fill mesh.

    Conclusion: the rung/diag Warren-truss bracing alone already prevents
    this artifact at current amplitude/material settings; enabling beams
    would currently be a no-op. Worth re-running the same rest-vs-late
    turning-angle comparison if AMPLITUDE or E_material change
    substantially, since that could plausibly push deformation into a
    regime where zig-zag actually starts to appear.

build_Tail_Ribbon_Connections (interior 'fill' springs):
    Each ribbon's actual 2D area (not just its two rails) is now filled
    with Halton-sampled interior points (_mesh_utils.fill_polygon_with_halton)
    and Delaunay-triangulated (_mesh_utils.triangulate_with_holes) into
    'fill' springs, controlled by case_params.json's geom.ribbon_fill_*
    keys. Their stiffness is NOT the rail/rung/diag flange-web formula used
    elsewhere in this file -- it's a separate tributary-area lattice-spring
    approximation:

        k_edge = ribbon_fill_stiffness_scale * E_material * extrude_depth
                 * A_edge / L_edge**2

    where A_edge is 1/3 of the area of each triangle sharing that edge
    (equal-split lumping), derived by equating discrete spring energy to
    continuum strain energy over the edge's tributary volume
    (A_edge * extrude_depth).

    Known limitations, not yet resolved:
      - A uniform triangular central-force lattice has a fixed effective
        Poisson's ratio around 1/3 and cannot independently target the
        near-incompressible nu=0.499 assumed by shear_Stiffness_From_Material
        for the rung/diag springs -- the ribbon is now a hybrid of two
        different stiffness derivations at its rail boundary vs. its
        interior fill, not one consistent continuum model.
      - Only genuinely new edges (touching an interior point, or a
        boundary-boundary diagonal not already covered by the existing
        outline/ray/rung/diag logic) get this formula; existing rail
        edges are untouched.

    Recommended calibration before trusting comparisons across
    ribbon_fill_density_multiplier: same approach as the beam calibration
    above -- isolate a single ribbon (no fluid, no other rays), clamp one
    end via target points, apply a known tip force, and compare deflection
    to the analytical cantilever formula. Sweep ribbon_fill_stiffness_scale
    until simulated and analytical deflection agree. See
    verify_ribbon_fill.py for a fluid-free geometry/mesh sanity check
    (degenerate-spring / runaway-count check across several density
    multipliers) -- it does not do this force calibration.

    STABILITY WARNING (confirmed, not hypothetical): at this example's
    defaults (dt=1e-3, E_material=0.74e6, ribbon_fill_density_multiplier=1.0),
    ribbon_fill_stiffness_scale=1.0 (the formula's literal, uncalibrated
    value) makes the simulation go fully NaN partway through the run --
    reproduced by running to Tfinal=5.0 and inspecting the .vtk point data:

        scale=1.0, 0.5:  diverges to NaN around t=3.3-3.6s
        scale=0.1, 0.01: clean (no NaN) through the full Tfinal=5.0

    BASE_CASE in run_case.py now defaults to scale=0.1 because of this.
    The likely mechanism: some 'fill' edges reach k on the order of the
    model's previously-stiffest element (k_Attach), and adding dozens of
    comparably-stiff edges pushes the explicit time-stepping scheme (fixed
    dt, no dt/k coupling) past its stability margin -- consistent with a
    slowly-growing resonance rather than an immediate blow-up (it takes
    several oscillation periods of the FREQUENCY=1.0 Hz actuation to
    appear). This 0.1 default is an empirically-found *starting point*,
    not a substitute for the tip-deflection calibration above, and it will
    need re-checking (rerun a case to Tfinal and grep its last .vtk frame
    for "nan") whenever E_material, dt, or ribbon_fill_density_multiplier
    change materially from these defaults.

    UPDATE after the ray-to-outline weld fix (build_Tail_Beams note above):
    re-ran the same check (Tfinal=5.0, ribbon_fill_density_multiplier=1.0)
    with the welded-junction geometry and scale=1.0 no longer diverges --
    clean through the full run, where it previously went NaN by t~3.3-3.6s.
    Consistent with the likely mechanism above: k_Attach (50*k_Spring, the
    old model's stiffest element by far) is gone entirely now, so there's
    no longer an isolated very-stiff spring for the explicit integrator to
    struggle with at that joint. NOT yet re-tested at higher
    ribbon_fill_density_multiplier (8.0+) -- BASE_CASE's scale=0.1 default
    is left as-is (still valid, just conservative) rather than bumped back
    to 1.0, since multiplier=1.0 stability doesn't guarantee it at higher
    multipliers. Worth re-sweeping scale at higher multipliers now that the
    joint itself is no longer the likely bottleneck.

ribbon_fill_density_multiplier -> actual point count (reference):
    n_points_for_this_ribbon = round(multiplier * ribbon_area / ds**2),
    computed PER RIBBON (outline + each ray) via
    _mesh_utils.fill_polygon_with_halton, then summed -- NOT a global
    point count set directly, and NOT the same as "multiplier x points
    along the boundary" (boundary point count is controlled entirely by
    _n_points_for_length/ds and does not change with this multiplier at
    all). At this example's defaults (ds=15.625mm, wall_thickness=3.96mm),
    the ribbons are thin relative to ds, so even multiplier=1 ("~ds
    density") gives very few points:
        multiplier= 1 ->   6 points total  (outline=5, ray0=1, rays1-3=0 --
                                             their area rounds down to 0)
        multiplier= 4 ->  27
        multiplier= 8 ->  56
        multiplier=16 -> 109
        multiplier=32 -> 220
        multiplier=64 -> 439
    Per-ribbon areas at these defaults: outline 1281.9mm^2, ray0 (widest)
    147.6mm^2, ray1 110.8mm^2, ray2 82.4mm^2, ray3 (thinnest) 49.7mm^2.
    Recompute via ribbonPolys + _mesh_utils._shoelace_area if
    wall_thickness/L/W/rayPositions change.

ribbon_fill_stiffness_scale -- what it represents (reference):
    The fill mesh is a LATTICE-SPRING approximation of a continuum, not a
    continuum-FEM one: k_edge = scale * E_material * extrude_depth *
    A_edge / L_edge**2 is a first-principles-but-approximate conversion
    from "microscopic spring stiffness" to "macroscopic effective
    modulus" (see the formula derivation in the
    build_Tail_Ribbon_Connections note above). `scale` corrects for that
    conversion's inherent approximation -- it is NOT a correction to
    E_material itself, and scale=1.0 has no special claim to being
    "correct," just "the formula's literal value."

    HARD LIMIT, independent of scale: a uniform lattice of purely axial
    (central-force) springs has an effective Poisson's ratio FIXED at
    ~1/3, regardless of spring stiffness (a classical lattice-elasticity
    result, the "Cauchy relation" for 2D central-force lattices -- this is
    a topology property, not a tuning parameter). Real silicone is
    near-incompressible, nu~0.499. So this mesh can NEVER reproduce real
    silicone's Poisson's ratio no matter how scale is calibrated -- scale
    can only ever approximately correct the effective Young's modulus.
    Worth knowing before chasing a "perfect" scale value: there is a
    ceiling on how well this specific modeling choice can match silicone,
    independent of calibration effort.

    Quantitative calibration procedure recommended, NOT YET BUILT: isolate
    a single ribbon (no fluid, no other rays), clamp one end rigidly
    (stiff target points), apply a known force at the free tip, and
    compare the simulated deflection to TIMOSHENKO beam theory (not plain
    Euler-Bernoulli/bending-only -- the fill mesh's actual physical role
    here is closer to shear resistance than bending, since bending is
    already dominated by the rail/flange springs, so a bending-only
    comparison likely won't be sensitive to scale at all):
        delta = F*L^3/(3*E*I) + F*L/(kappa*G*A)
    with kappa~5/6 (rectangular cross-section shear correction) and
    G = E/(2*(1+nu)) using silicone's REAL nu~0.499 (not the mesh's own
    ~1/3). Solve for the scale value that makes simulated and analytical
    deflection agree. Mirrors NACA0020_Wing/verify_structure.py's existing
    fluid-free beam calibration pattern (Newton solve on IB2d's own force
    routines, no need to run the actual fluid solver).

Incomplete-rename / dropped-fallback bugs (FOUND + FIXED):
    Found while investigating the beams question above --
    _Build_Tail_Geometry_Ribbon raised immediately on any call:
        TypeError: fill_polygon_with_halton() got an unexpected keyword
        argument 'inner_loops'
    Two separate issues, both leftovers from in-progress edits that
    weren't fully propagated:

    1. ribbonPolys' dict keys (exterior->exterior_loop, holes->
       inner_loops) and _mesh_utils.triangulate_with_holes's *positional*
       call had been renamed consistently throughout FinRay_Geom.py, but
       fill_polygon_with_halton is called with inner_loops= as a KEYWORD
       argument (lines ~827, 836), and _mesh_utils.py's own parameter was
       still named `holes` -- hence the TypeError.
       triangulate_with_holes itself wasn't broken (called positionally,
       so its internal parameter name doesn't matter to callers), but was
       left inconsistent with fill_polygon_with_halton. Fixed by
       completing the rename in _mesh_utils.py: both functions' parameter
       name, plus their internal variable names/docstrings
       (_far_from_boundary's `holes` param too).

    2. _Build_Tail_Geometry_Ribbon's own
       `fill_spacing = ds if fill_spacing is None else fill_spacing`
       normalization line had been deleted -- apparently on the
       assumption that resolving case_params.json's default earlier, in
       FinRay_Geom() itself (`fill_spacing = _geom.get('ribbon_fill_spacing',
       ds)`), made it redundant. But the function's own docstring still
       promises "fill_spacing: ... None falls back to ds", and at least
       one other caller, give_Me_Immersed_Boundary_Geometry, still calls
       it with `fill_spacing=_geom.get('ribbon_fill_spacing', None)` --
       i.e. still relies on that internal fallback, and would pass None
       straight through to fill_polygon_with_halton (-> ValueError:
       "needs spacing or n_target") whenever ribbon_fill_spacing isn't
       set in case_params.json (the common case, since it defaults to
       null there). Restored the line.

    Verified via a fresh FinRay_Geom() run + verify_ribbon_fill.py:
    identical point/spring counts to before these two bugs were
    introduced (298 total points for the mult32.00_scale1.00 config,
    matching that case's already-archived results exactly), confirming
    the fix is behavior-preserving, not a new behavior change.