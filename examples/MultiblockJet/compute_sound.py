"""
Compute microphone pressure histories and SPL from the baseline
Mach 1.3 MultiblockJet FWH probe data.

Formulation: kim_thesis.pdf sections 4.2, 6.2; thesis.pdf section 8.3.
Reuses postprocess.compute_sound (FWH integration), fwhsolver.windowed_fft
(Blackman-windowed FFT, 50% overlap, returns freq and RMS-averaged p_hat),
and postprocess.get_SPL (converts p_hat to SPL/OASPL in dB).
"""
import argparse
import os
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from magudi_utils import plot3dnasa as p3d
from magudi_utils import fwhsolver as fwh
from postprocess import (compute_sound, get_SPL,
                         polynomial_pressure_mean, moving_average)


# bc.dat: fwh1..fwh4 PROBE patches sit at these radial indices in blocks 2-5.
SURFACE_TO_RADIAL_INDEX = {
    'fwh1': 155,
    'fwh2': 172,
    'fwh3': 183,
    'fwh4': 192,
}


AUTHOR_STYLES = {
    'Kim':    dict(color='magenta', ls='-.', lw=1.3),
    'Samimy': dict(color='black',   ls='--', lw=1.3),
    'Ram':    dict(color='blue',    ls=(0, (1, 1)), lw=1.3),
}


def fwh_surface_radius(grid_file, surface):
    """Return the physical radius of the FWH surface `surface` (fwh1..fwh4),
    read from block 1 (0-indexed) of `grid_file` at the mapped radial index."""
    idx = SURFACE_TO_RADIAL_INDEX[surface]
    g = p3d.Grid(grid_file).set_subzone(1, [idx, 0, 0], [idx, -1, 0]).load()
    r = np.mean(np.sqrt(g.xyz[0][0, :, 0, 0] ** 2 +
                        g.xyz[0][0, :, 0, 1] ** 2))
    return float(r)


def load_mike_pressures(dir, probe_name, num_mikes):
    """Read the mike_{probe_name}_##.dat files back into p[nsteps, num_mikes]."""
    columns = []
    for i in range(num_mikes):
        data = np.loadtxt('%s/mike_%s_%02d.dat' % (dir, probe_name, i + 1))
        columns.append(data[:, 1])
    return np.array(columns).T


def _overlay_literature_spl(ax, literature_h5, station):
    with h5py.File(literature_h5, 'r') as f:
        for author in f:
            key = '%s/SPL/%s' % (author, station)
            if key not in f:
                continue
            arr = f[key][...]
            style = AUTHOR_STYLES.get(author, dict(lw=1.))
            ax.plot(arr[:, 0], arr[:, 1], label=author, **style)


def diagnose_stationarity(p_raw, p, p_raw_avg, p_avg,
                          num_windows, dt, mach_number,
                          distance, theta, literature_h5, out_path):
    """Split the (detrended) mike record in half, compute SPL on each half, and
    overlay with the full-signal SPL, the full raw-signal SPL, the SPL of the
    mike-averaged pressure (raw and detrended), and literature. If the two
    halves overlay, the signal is statistically stationary; if the first half
    sits higher (especially at low St), an initial transient is contaminating
    the record. The raw-vs-detrended comparison isolates the effect of
    polynomial mean-drift subtraction on the low-St band.

    `p_raw_avg` / `p_avg` are the mike-averaged raw and detrended pressures
    (1-D length nsteps), where the detrending is a polynomial fit to the
    averaged signal itself (not the average of the per-mike detrended signals).
    windowed_fft needs 2D input, so they are reshaped internally."""
    n = p.shape[0]
    kw_fft = dict(num_windows=num_windows, dt=dt)
    kw_spl = dict(dt=dt, mach_number=mach_number, distance=distance)
    freq_raw,     phat_raw     = fwh.windowed_fft(p_raw,     **kw_fft)
    freq_full,    phat_full    = fwh.windowed_fft(p,         **kw_fft)
    freq_raw_avg, phat_raw_avg = fwh.windowed_fft(p_raw_avg[:, None], **kw_fft)
    freq_avg,     phat_avg     = fwh.windowed_fft(p_avg[:, None],     **kw_fft)
    freq1,        phat1        = fwh.windowed_fft(p[:n // 2], **kw_fft)
    freq2,        phat2        = fwh.windowed_fft(p[n // 2:], **kw_fft)
    St_raw,     SPL_raw,     _ = get_SPL(freq_raw,     phat_raw,     **kw_spl)
    St_full,    SPL_full,    _ = get_SPL(freq_full,    phat_full,    **kw_spl)
    St_raw_avg, SPL_raw_avg, _ = get_SPL(freq_raw_avg, phat_raw_avg, **kw_spl)
    St_avg,     SPL_avg,     _ = get_SPL(freq_avg,     phat_avg,     **kw_spl)
    St1,        SPL1,        _ = get_SPL(freq1,        phat1,        **kw_spl)
    St2,        SPL2,        _ = get_SPL(freq2,        phat2,        **kw_spl)

    station = '%gD%gdeg' % (distance, theta)
    fig, ax = plt.subplots(1, 1, figsize=(6, 4.5))
    ax.semilogx(St_raw,     SPL_raw,     color='gray',   ls=':',  lw=1.6,
                label='Full (raw)')
    ax.semilogx(St_full,    SPL_full,    color='red',    ls='-',  lw=1.6,
                label='Full (detrended)')
    ax.semilogx(St_raw_avg, SPL_raw_avg, color='purple', ls=':',  lw=1.3,
                label='Mike avg (raw)')
    ax.semilogx(St_avg,     SPL_avg,     color='purple', ls='-',  lw=1.3,
                label='Mike avg (detrended)')
    ax.semilogx(St1, SPL1, color='orange', ls='--', lw=1.3, label='1st half')
    ax.semilogx(St2, SPL2, color='green',  ls='--', lw=1.3, label='2nd half')
    _overlay_literature_spl(ax, literature_h5, station)
    ax.set_xlabel(r'$St_D$')
    ax.set_ylabel('SPL (dB)')
    ax.set_title(r'Stationarity diagnostic — $d=%gD$, $\theta=%g^\circ$' %
                 (distance, theta))
    # ax.set_xlim(0.05, 5.0)
    ax.set_xlim(0.005, 5.0)
    ax.set_ylim(0.0, 110.0)
    ax.set_yticks(np.arange(0, 111, 10))
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(loc='best', frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_mike_history(p, dt, distance, theta, out_path,
                      poly_mean=None, moving_average=None,
                      p_avg=None, poly_mean_avg=None,
                      moving_average_avg=None):
    """Plot per-mike pressure time history as (num_mikes, 1) subplots.

    If `poly_mean` (blue dashed) and/or `moving_average` (red solid) are
    supplied, overlay them on each subplot for slow-drift diagnosis.

    If `p_avg` is supplied, append one extra subplot at the bottom showing the
    mike-averaged pressure (black), with optional `poly_mean_avg` (blue dashed)
    and `moving_average_avg` (red solid) overlays."""
    nsteps, num_mikes = p.shape
    t = np.arange(nsteps) * dt
    has_avg = p_avg is not None
    n_rows = num_mikes + (1 if has_avg else 0)
    fig, axes = plt.subplots(n_rows, 1, sharex=True,
                             figsize=(10, 1.2 * n_rows))
    if n_rows == 1:
        axes = [axes]
    for i in range(num_mikes):
        ax = axes[i]
        ax.plot(t, p[:, i], 'k-', lw=0.6)
        if poly_mean is not None:
            ax.plot(t, poly_mean[:, i], 'b--', lw=1.0)
        if moving_average is not None:
            ax.plot(t, moving_average[:, i], 'r-', lw=1.0)
        ax.set_ylabel('mike %d' % (i + 1), fontsize=8)
        ax.tick_params(labelsize=8)
    if has_avg:
        ax = axes[-1]
        ax.plot(t, p_avg, 'k-', lw=0.6)
        if poly_mean_avg is not None:
            ax.plot(t, poly_mean_avg, 'b--', lw=1.0)
        if moving_average_avg is not None:
            ax.plot(t, moving_average_avg, 'r-', lw=1.0)
        ax.set_ylabel('mike avg', fontsize=8)
        ax.tick_params(labelsize=8)
    axes[-1].set_xlabel(r'$t\,a_\infty / D$')
    axes[0].set_title(r'Mike histories — $d=%gD$, $\theta=%g^\circ$' %
                      (distance, theta))
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_SPL(spl_file, distance, theta, literature_h5, out_path):
    """Plot SPL(St) from `spl_file` against literature data from `literature_h5`
    at the station keyed by (distance, theta), e.g. (94, 30) -> '94D30deg'."""
    station = '%gD%gdeg' % (distance, theta)

    fig, ax = plt.subplots(1, 1, figsize=(6, 4.5))
    d = np.loadtxt(spl_file)
    ax.semilogx(d[:, 0], d[:, 1], 'r-', lw=1.6, label='Current')
    _overlay_literature_spl(ax, literature_h5, station)
    ax.set_xlabel(r'$St_D$')
    ax.set_ylabel('SPL (dB)')
    ax.set_title(r'$d=%gD$, $\theta=%g^\circ$' % (distance, theta))
    ax.set_xlim(0.05, 5.0)
    # ax.set_xlim(0.005, 5.0)
    ax.set_ylim(0.0, 110.0)
    ax.set_yticks(np.arange(0, 111, 10))
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(loc='best', frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prefix', default='MultiblockJet')
    p.add_argument('--surface', default='fwh1',
                   choices=sorted(SURFACE_TO_RADIAL_INDEX.keys()))
    p.add_argument('--distance', type=float, required=True,
                   help='slant distance to microphone ring')
    p.add_argument('--theta', type=float, required=True,
                   help='polar angle in degrees')
    p.add_argument('--num-mikes', type=int, default=8)
    p.add_argument('--x0', type=float, nargs=3, default=[0.0, 0.0, -1.8],
                   metavar=('X', 'Y', 'Z'),
                   help='jet origin (default reflects nozzle-exit offset)')
    p.add_argument('--probe-dt', type=float, default=35 * 1.2e-3,
                   help='nondim time between probe samples')
    p.add_argument('--mach', type=float, default=1.3)
    p.add_argument('--num-windows', type=int, default=5)
    p.add_argument('--poly-degree', type=int, default=3,
                   help='degree of least-squares polynomial fit to per-mike '
                        'pressure history; subtracted from p before FFT to '
                        'remove slow mean drift')
    p.add_argument('--ma-frac', type=float, default=0.05,
                   help='moving-average window as a fraction of the record '
                        'length (used only for the mike-history overlay)')
    p.add_argument('--load-mike', action='store_true',
                   help='skip the FWH integration and load existing '
                        'mike_<surface>_##.dat files instead; still re-runs '
                        'windowed_fft and re-writes spl/oaspl outputs')
    p.add_argument('--diagnose', action='store_true',
                   help='stationarity diagnostic: split mike record in half, '
                        'plot SPL of each half + full signal against literature. '
                        'Assumes mike_<surface>_##.dat already exist in '
                        '--spl-dir; skips FWH, spl/oaspl outputs, and plot_SPL.')
    p.add_argument('--spl-dir', default='.',
                   help='directory for all outputs written by this script: '
                        'mike_<surface>_##.dat, spl_<surface>.dat, '
                        'oaspl_<surface>.dat, and (when --figure is not set) '
                        'the SPL/diagnostic PNGs')
    p.add_argument('--literature', default='literature.h5',
                   help='HDF5 file with literature SPL data '
                        '(layout: /<Author>/SPL/{94D30deg,44D90deg})')
    p.add_argument('--figure', default=None,
                   help='output path for the SPL figure (default: '
                        'spl_<surface>_<distance>D_<theta>deg.png)')
    return p.parse_args()


def main():
    args = parse_args()

    os.makedirs(args.spl_dir, exist_ok=True)
    spl_file = os.path.join(args.spl_dir, 'spl_%s.dat' % args.surface)
    oaspl_file = os.path.join(args.spl_dir, 'oaspl_%s.dat' % args.surface)

    if not args.load_mike:
        grid_file = '%s.xyz' % args.prefix

        probe_r = fwh_surface_radius(grid_file, args.surface)
        print('%s: radial index %d -> r/D = %.4f' %
              (args.surface, SURFACE_TO_RADIAL_INDEX[args.surface], probe_r))

        compute_sound(args.prefix, args.x0, args.probe_dt, args.distance,
                      args.theta, args.surface, probe_r, out_dir=args.spl_dir)

    p_raw = load_mike_pressures(args.spl_dir, args.surface, args.num_mikes)

    poly_mean = polynomial_pressure_mean(p_raw, degree=args.poly_degree)
    ma_window = max(1, int(round(args.ma_frac * p_raw.shape[0])))
    ma_mean = moving_average(p_raw, ma_window)

    # Mike-averaged pressure and its own polynomial / moving-average fits
    # (fit on the averaged signal, not the average of the per-mike fits).
    p_raw_avg = p_raw.mean(axis=1)
    poly_mean_avg = polynomial_pressure_mean(
        p_raw_avg[:, None], degree=args.poly_degree)[:, 0]
    ma_mean_avg = moving_average(p_raw_avg[:, None], ma_window)[:, 0]

    history_fig = os.path.join(
        args.spl_dir, 'mike_history_%s_%gD_%gdeg.png' %
        (args.surface, args.distance, args.theta))
    plot_mike_history(p_raw, args.probe_dt, args.distance, args.theta,
                      history_fig,
                      poly_mean=poly_mean, moving_average=ma_mean,
                      p_avg=p_raw_avg, poly_mean_avg=poly_mean_avg,
                      moving_average_avg=ma_mean_avg)
    print('wrote %s' % history_fig)

    # p = p_raw - poly_mean
    p = p_raw - ma_mean
    p_avg = p_raw_avg - ma_mean_avg

    if args.diagnose:
        diag_fig = args.figure or os.path.join(
            args.spl_dir, 'diagnose_%s_%gD_%gdeg.png' %
            (args.surface, args.distance, args.theta))
        diagnose_stationarity(p_raw, p, p_raw_avg, p_avg,
                              args.num_windows, args.probe_dt,
                              args.mach, args.distance, args.theta,
                              args.literature, diag_fig)
        print('wrote %s' % diag_fig)
        return

    freq, p_hat = fwh.windowed_fft(p, num_windows=args.num_windows,
                                   dt=args.probe_dt)
    St, SPL, OASPL = get_SPL(freq, p_hat, dt=args.probe_dt,
                             mach_number=args.mach, distance=args.distance)

    np.savetxt(spl_file, np.column_stack([St, SPL]), fmt='%+.18E')
    with open(oaspl_file, 'a') as f:
        f.write('%+.6E %+.6E\n' % (args.theta, OASPL))
    print('%s @ theta=%g deg, d=%g: OASPL = %.3f dB' %
          (args.surface, args.theta, args.distance, OASPL))

    figure = args.figure or os.path.join(
        args.spl_dir, 'spl_%s_%gD_%gdeg.png' %
        (args.surface, args.distance, args.theta))
    plot_SPL(spl_file, args.distance, args.theta, args.literature, figure)
    print('wrote %s' % figure)


if __name__ == '__main__':
    main()
