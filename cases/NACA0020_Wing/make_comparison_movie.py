'''
Side-by-side vorticity movie: a SoRoLily run (its own NACA0020.mp4, made by
gvs_tests.jl) next to an IB2d run (its viz_IB2d/Omega.*.vtk), frame by frame
at the same instants, same colour scale, same zoom.

    python make_comparison_movie.py SOROLILY_CASE_DIR IB2D_RUN_DIR OUT.mp4 \
        "left label" "right label" [--snapshots t1,t2,...]

Both codes save a frame every 0.05 s (SoRoLily: fps = 20 per simulated second,
from t = 0; IB2d: field_interval = 0.05 s), so frame k is t = 0.05 k in both.
SoRoLily colours curl(u) * L_c / U_c, i.e. the vorticity times c/U = 1 s, so
it is numerically the vorticity in 1/s, which is what IB2d writes; both are
drawn on [-10, 10] with ColorBrewer RdBu (red negative, blue positive), as
SoRoLily's Plots :RdBu_11. SoRoLily's frame is cropped to its plotted domain
(found from the frame itself) and placed on the same 1.2 x 1.2 m axes.
With --snapshots, also writes OUT_snapshots.png at those times (seconds).
Needs ffmpeg.
'''

import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE / 'reference'))
import sorolily as sl
import wing_case
from postprocess import read_scalar_vtk

FRAME_DT = 0.05
XLIM, YLIM = (0.25, 1.2), (0.35, 0.85)     # m: the wing and its near wake


def soro_frames(case_dir, tmp):
    subprocess.run(['ffmpeg', '-v', 'error', '-y', '-i', str(Path(case_dir) / 'NACA0020.mp4'),
                    str(Path(tmp) / 'soro_%04d.png')], check=True)
    return sorted(Path(tmp).glob('soro_*.png'))


def domain_box(img):
    '''Pixel box of SoRoLily's plotted domain: the non-white area left of the colour bar.'''
    rgb = img[..., :3]
    if rgb.max() > 1.5:
        rgb = rgb / 255.0
    w = rgb.shape[1]
    non_white = np.any(rgb[:, : int(0.8 * w)] < 0.985, axis=2)
    rows, cols = np.where(non_white)
    return rows.min(), rows.max(), cols.min(), cols.max()


def ib2d_frame(run_dir, k, hist, n_out):
    f = Path(run_dir) / 'viz_IB2d' / f'Omega.{k:04d}.vtk'
    if not f.exists():
        return None, None
    w = read_scalar_vtk(f)
    i = int(np.argmin(np.abs(hist['t'] - k * FRAME_DT)))
    return w, hist['pos'][i, :n_out].astype(float)


def draw(ax_l, ax_r, soro_img, box, w, outline, Lx, labels, t, cmap, levels):
    r0, r1, c0, c1 = box
    ax_l.imshow(soro_img[r0:r1 + 1, c0:c1 + 1], extent=(0, sl.DOMAIN[0], 0, sl.DOMAIN[1]))
    ny, nx = w.shape
    xs, ys = np.arange(nx) * Lx / nx, np.arange(ny) * Lx / ny
    cf = ax_r.contourf(xs, ys, np.clip(w, -10, 10), levels=levels, cmap=cmap, extend='both')
    o = np.vstack([outline, outline[:1]])
    ax_r.fill(o[:, 0], o[:, 1], color='white', zorder=3)
    ax_r.plot(o[:, 0], o[:, 1], 'k-', lw=1.5, zorder=4)
    for ax, lab in ((ax_l, labels[0]), (ax_r, labels[1])):
        ax.set_xlim(*XLIM); ax.set_ylim(*YLIM); ax.set_aspect('equal')
        ax.set_title(lab, fontsize=12); ax.set_xlabel('x [m]')
    ax_l.set_ylabel('y [m]')
    return cf


def main(argv):
    snaps = None
    if '--snapshots' in argv:
        i = argv.index('--snapshots'); snaps = [float(v) for v in argv[i + 1].split(',')]; del argv[i:i + 2]
    soro_dir, run_dir, out, lab_l, lab_r = argv[:5]
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.image import imread
    hist = np.load(Path(run_dir) / 'viz_IB2d' / 'wing_history.npz')
    n_out = int(np.load(Path(run_dir) / 'wing_mesh.npz')['n_outline'])
    from run_records import read_input2d
    Lx = read_input2d(Path(run_dir) / 'input2d')['Lx']
    flow = __import__('json').loads((Path(run_dir) / 'case_params.json').read_text())['flow']
    cmap, levels = plt.get_cmap('RdBu'), np.linspace(-10, 10, 41)
    with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR')) as tmp:
        frames = soro_frames(soro_dir, tmp)
        box = domain_box(imread(frames[len(frames) // 2]))
        n = len(frames)
        made = 0
        for k in range(n):
            w, outline = ib2d_frame(run_dir, k, hist, n_out)
            if w is None or not np.all(np.isfinite(w)):
                break
            t = k * FRAME_DT
            fig, (al, ar) = plt.subplots(1, 2, figsize=(16, 5.6))
            cf = draw(al, ar, imread(frames[k]), box, w, outline, Lx, (lab_l, lab_r), t, cmap, levels)
            fig.colorbar(cf, ax=[al, ar], shrink=0.85, label='vorticity [1/s] (= ω c/U)', ticks=[-10, -5, 0, 5, 10])
            aoa = np.degrees(wing_case.aoa_rad(t, flow))
            fig.suptitle(f't = {t:5.2f} s     angle of attack {aoa:4.1f}°', fontsize=13)
            fig.savefig(Path(tmp) / f'cmp_{k:04d}.png', dpi=90); plt.close(fig)
            made += 1
        subprocess.run(['ffmpeg', '-v', 'error', '-y', '-framerate', '20', '-i', str(Path(tmp) / 'cmp_%04d.png'),
                        '-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2', '-c:v', 'libx264', '-pix_fmt', 'yuv420p', out],
                       check=True)
        print(f'wrote {out}: {made} frames ({made * FRAME_DT:.2f} s)')
        if snaps:
            fig, axs = plt.subplots(len(snaps), 2, figsize=(14, 3.6 * len(snaps)))
            for row, ts in zip(np.atleast_2d(axs), snaps):
                k = int(round(ts / FRAME_DT))
                w, outline = ib2d_frame(run_dir, k, hist, n_out)
                cf = draw(row[0], row[1], imread(frames[k]), box, w, outline, Lx,
                          (f'{lab_l}, t = {ts:g} s', f'{lab_r}, t = {ts:g} s'), ts, cmap, levels)
            fig.colorbar(cf, ax=axs, shrink=0.6, label='vorticity [1/s] (= ω c/U)', ticks=[-10, -5, 0, 5, 10])
            png = os.path.splitext(out)[0] + '_snapshots.png'
            fig.savefig(png, dpi=110); plt.close(fig)
            print('wrote', png)


if __name__ == '__main__':
    if len(sys.argv) < 6:
        sys.exit(__doc__)
    main(sys.argv[1:])
