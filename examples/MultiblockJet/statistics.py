"""
Compute Favre-averaged streamwise velocity on the centerline (r/D = 0)
and nozzle lipline (r/D = 0.5) from <input_dir>/<prefix>.mean.q, write
a three-column text file, and produce a 2x2 comparison figure:

                mean                       rms
    +-------------------------+ +-------------------------+
    | Centerline: sim + lit   | | Centerline: literature  |
    +-------------------------+ +-------------------------+
    | Lipline:    sim + lit   | | Lipline:    literature  |
    +-------------------------+ +-------------------------+

The literature curves are loaded from an HDF5 file with the layout
    /<Author>/<station>/{mean, rms}     (dataset)
    /<Author>/@reference                (attribute)
as produced by build_literature_h5.py.
"""
import argparse
import os
import warnings
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from magudi_utils import plot3dnasa as p3d
from magudi_utils.RoundJet import getCenterlineSubarray, mergeMeanSolutions, \
    getCenterlineRMSFluctuations


# --- copied verbatim from examples/MultiblockJet/postprocess.py ---
def extract_const_r(g, f, r=0.5):
    n = f.get_size()
    g.set_subzone(1, [0, 0, 0], [-1, -1, 0]).load()
    idx = np.argmin(np.abs(np.mean(np.sqrt(g.xyz[0][:,:,0,0] ** 2 +
                                           g.xyz[0][:,:,0,1] ** 2),
                                   axis=1) - r))
    args = dict()
    if type(f) == p3d.Function:
        args.update(ncomponents=f.ncomponents)
    fe = type(f)(**args).set_size([n[0][2], 4 * (n[1][1] - 1) + 1, 1], True)
    for i in range(1, 5):
        f.set_subzone(i, [idx, 0, 0], [idx, -2, -1]).load()
        for j in range(f[0].shape[-1]):
            fe[0][:,(i-1)*(n[1][1]-1):i*(n[1][1]-1),0,j] = f[0][0,:,:,j].T
        fe[0][:,-1,:,:] = fe[0][:,0,:,:]
    return fe
# ------------------------------------------------------------------


def merge_means(prefix, mean_files):
    out = '%s.mean.q' % prefix
    if os.path.exists(out):
        raise RuntimeError('%s already exists' % out)
    mergeMeanSolutions(prefix, mean_files).save(out)


def centerline_favre_u3(prefix):
    z, Q = getCenterlineSubarray('%s.xyz' % prefix, '%s.mean.q' % prefix)
    return z, Q[:, 3] / Q[:, 0]


def lipline_radius(g, r=0.5, r_tol=1e-10):
    g.set_subzone(1, [0, 0, 0], [-1, -1, 0]).load()
    rr = np.mean(np.sqrt(g.xyz[0][:, :, 0, 0] ** 2 +
                         g.xyz[0][:, :, 0, 1] ** 2), axis=1)
    r_hit = rr[np.argmin(np.abs(rr - r))]
    print('lipline: nearest index r/D = %.4f (target %.4f)' % (r_hit, r))
    if abs(r_hit - r) > r_tol:
        warnings.warn('lipline nearest index r/D = %.4f deviates from target '
                      '%.4f by more than %g' % (r_hit, r, r_tol))
    return r_hit


def lipline_favre_u3(prefix, r=0.5):
    g = p3d.Grid('%s.xyz' % prefix)
    q = p3d.fromfile('%s.mean.q' % prefix)
    fe = extract_const_r(g, q, r=r)
    ring = fe[0][:, :-1, 0, :]
    rhobar     = np.mean(ring[:, :, 0], axis=1)
    rho_u3_bar = np.mean(ring[:, :, 3], axis=1)
    z = g.set_subzone(0, [0, 0, 0], [0, 0, -1]).load().xyz[0][0, 0, :, 2]
    return z, rho_u3_bar / rhobar


def lipline_rms_u3(prefix, soln_files, r=0.5, gamma=1.4):
    g = p3d.Grid('%s.xyz' % prefix)
    ring_m = extract_const_r(g, p3d.fromfile('%s.mean.q' % prefix), r=r)[0]
    u3m = ring_m[:, :-1, 0, 3] / ring_m[:, :-1, 0, 0]
    acc = np.zeros_like(u3m)
    for f in soln_files:
        ring = extract_const_r(g, p3d.Solution(f), r=r)[0]
        u3 = ring[:, :-1, 0, 3] / ring[:, :-1, 0, 0]
        acc += (u3 - u3m) ** 2
    rms = np.sqrt(np.mean(acc / len(soln_files), axis=1))
    z = g.set_subzone(0, [0, 0, 0], [0, 0, -1]).load().xyz[0][0, 0, :, 2]
    return z, rms


AUTHOR_STYLES = {
    'Kim':     dict(color='magenta', ls=(0, (6, 2)), lw=1.3),
    'Bodony':  dict(color='green',   ls='--',        lw=1.3),
    'Mendez':  dict(color='blue',    ls='-.',        lw=1.3),
    'Bridges': dict(color='cyan',    marker='o', mfc='none', ls='none', ms=5),
    'Samimy':  dict(color='black',   marker='s', mfc='none', ls='none', ms=5),
}


def overlay_literature(ax, h5path, station, statistic):
    if h5path is None or not os.path.exists(h5path):
        return
    with h5py.File(h5path, 'r') as f:
        for author in f:
            g_a = f[author]
            if station not in g_a or statistic not in g_a[station]:
                continue
            d = g_a[station][statistic][...]
            style = AUTHOR_STYLES.get(author, dict(lw=1.))
            ax.plot(d[:, 0], d[:, 1], label=author, **style)


XLIM = {
    ('centerline', 'mean'): (1.0, 24.0),   # Fig 8.5(a)
    ('centerline', 'rms'):  (1.5, 15.0),   # Fig 8.5(b)
    ('lipline',    'mean'): (0.0, 25.0),   # Fig 8.6
    ('lipline',    'rms'):  (0.0, 25.0),   # user-specified
}


def make_figure(x_sim, u3_centerline, u3_lipline, literature_h5, out_path,
                shift=0., rms_centerline=None, rms_lipline=None):
    current_label = r'Current ($x_s=%g$)' % shift
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))

    # centerline mean : simulation + literature
    ax = axes[0, 0]
    ax.plot(x_sim, u3_centerline, 'r-', lw=1.6, label=current_label)
    overlay_literature(ax, literature_h5, 'centerline', 'mean')
    ax.set_ylabel(r'$\tilde{u}_3/U_j$')
    ax.set_title(r'Centerline ($r/D=0$) — mean')
    ax.set_xlabel(r'$(x_3 + x_s)/D$')
    ax.set_xlim(*XLIM[('centerline', 'mean')])

    # centerline rms : simulation (if available) + literature
    ax = axes[0, 1]
    if rms_centerline is not None:
        ax.plot(x_sim, rms_centerline, 'r-', lw=1.6, label=current_label)
    overlay_literature(ax, literature_h5, 'centerline', 'rms')
    ax.set_ylabel(r'$(\overline{u_3^{\prime\prime} u_3^{\prime\prime}})^{1/2}/U_j$')
    ax.set_title(r'Centerline ($r/D=0$) — rms')
    ax.set_xlabel(r'$(x_3 + x_s)/D$')
    ax.set_xlim(*XLIM[('centerline', 'rms')])

    # lipline mean : simulation + literature
    ax = axes[1, 0]
    ax.plot(x_sim, u3_lipline, 'r-', lw=1.6, label=current_label)
    overlay_literature(ax, literature_h5, 'lipline', 'mean')
    ax.set_ylabel(r'$\tilde{u}_3/U_j$')
    ax.set_title(r'Nozzle lipline ($r/D=0.5$) — mean')
    ax.set_xlabel(r'$(x_3 + x_s)/D$')
    ax.set_xlim(*XLIM[('lipline', 'mean')])

    # lipline rms : simulation (if available) + literature
    ax = axes[1, 1]
    if rms_lipline is not None:
        ax.plot(x_sim, rms_lipline, 'r-', lw=1.6, label=current_label)
    overlay_literature(ax, literature_h5, 'lipline', 'rms')
    ax.set_ylabel(r'$(\overline{u_3^{\prime\prime} u_3^{\prime\prime}})^{1/2}/U_j$')
    ax.set_title(r'Nozzle lipline ($r/D=0.5$) — rms')
    ax.set_xlabel(r'$(x_3 + x_s)/D$')
    ax.set_xlim(*XLIM[('lipline', 'rms')])

    for ax in axes.ravel():
        ax.grid(True, alpha=0.3)
        if ax.get_legend_handles_labels()[1]:
            ax.legend(loc='best', frameon=False, fontsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main(input_dir='.', output_dir='.', prefix='MultiblockJet', shift=0.,
         literature_h5=None, mach_number=1.3, gamma=1.4, r_lipline=0.5,
         mean_files=None, start=None, end=None, freq=35):
    in_prefix = os.path.join(input_dir, prefix)
    out_txt   = os.path.join(output_dir, '%s.favre_u3.txt' % prefix)
    out_png   = os.path.join(output_dir, '%s.stats.png'    % prefix)

    if mean_files is not None:
        merge_means(in_prefix, mean_files)

    z_c, u3_c = centerline_favre_u3(in_prefix)
    z_l, u3_l = lipline_favre_u3(in_prefix, r=r_lipline)
    assert np.allclose(z_c, z_l)

    T_ratio = 1. / (1. + 0.5 * (gamma - 1.) * mach_number ** 2)
    U_j = mach_number * np.sqrt(T_ratio)

    x        = z_c + shift
    u3_c_hat = u3_c / U_j
    u3_l_hat = u3_l / U_j

    rms_c_hat, rms_l_hat = None, None
    if start is not None and end is not None:
        lipline_radius(p3d.Grid('%s.xyz' % in_prefix), r=r_lipline)
        soln_files = [os.path.join(input_dir, '%s.%08d.q' % (prefix, t))
                      for t in range(start, end + 1, freq)]
        rms_c_hat = getCenterlineRMSFluctuations(
            '%s.mean.q' % in_prefix, soln_files,
            ratioOfSpecificHeats=gamma)[:, 3] / U_j
        z_lr, rms_l = lipline_rms_u3(in_prefix, soln_files, r=r_lipline,
                                     gamma=gamma)
        assert np.allclose(z_lr, z_c)
        rms_l_hat = rms_l / U_j

    if rms_c_hat is None:
        np.savetxt(
            out_txt, np.column_stack([x, u3_c_hat, u3_l_hat]),
            header='(x3+xs)/D    u_tilde_3/U_j (r/D=0)    u_tilde_3/U_j (r/D=0.5)',
            fmt='%+.15E')
    else:
        np.savetxt(
            out_txt, np.column_stack([x, u3_c_hat, u3_l_hat,
                                      rms_c_hat, rms_l_hat]),
            header='(x3+xs)/D    u_tilde_3/U_j (r/D=0)    '
                   'u_tilde_3/U_j (r/D=0.5)    rms_u3/U_j (r/D=0)    '
                   'rms_u3/U_j (r/D=0.5)',
            fmt='%+.15E')

    make_figure(x, u3_c_hat, u3_l_hat, literature_h5, out_png, shift=shift,
                rms_centerline=rms_c_hat, rms_lipline=rms_l_hat)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input-dir',  default='.',
                        help='directory holding <prefix>.xyz and <prefix>.mean.q')
    parser.add_argument('--output-dir', default='.',
                        help='directory for the output .txt and .png files')
    parser.add_argument('--prefix',     default='MultiblockJet',
                        help='file prefix used for both inputs and outputs')
    parser.add_argument('--shift',      type=float, default=0.,
                        help='shifting distance x_s/D applied to the axial coordinate '
                             '(x_s=1.8 reproduces the thesis figures)')
    parser.add_argument('--literature', default=None,
                        help='path to an HDF5 file with literature statistics '
                             '(layout: /<Author>/{centerline,lipline}/{mean,rms})')
    parser.add_argument('--mean-files', nargs='+', default=None,
                        help='list of mean .q files to arithmetic-average into '
                             '<prefix>.mean.q (fails if that file already exists)')
    parser.add_argument('--start', type=int, default=None,
                        help='first snapshot timestep for the rms computation')
    parser.add_argument('--end', type=int, default=None,
                        help='last snapshot timestep for the rms computation')
    parser.add_argument('--freq', type=int, default=35,
                        help='timestep stride between snapshots (default 35)')
    args = parser.parse_args()
    main(input_dir=args.input_dir, output_dir=args.output_dir,
         prefix=args.prefix, shift=args.shift, literature_h5=args.literature,
         mean_files=args.mean_files, start=args.start, end=args.end,
         freq=args.freq)
