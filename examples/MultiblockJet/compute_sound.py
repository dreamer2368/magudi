"""
Compute microphone pressure histories and SPL from the baseline
Mach 1.3 MultiblockJet FWH probe data.

Formulation: kim_thesis.pdf sections 4.2, 6.2; thesis.pdf section 8.3.
Reuses postprocess.compute_sound (FWH integration) and postprocess.windowed_fft
(Blackman-windowed FFT, 50% overlap, returns Strouhal, SPL, OASPL).
"""
import argparse
import numpy as np
from magudi_utils import plot3dnasa as p3d
from postprocess import compute_sound, windowed_fft


# bc.dat: fwh1..fwh4 PROBE patches sit at these radial indices in blocks 2-5.
SURFACE_TO_RADIAL_INDEX = {
    'fwh1': 155,
    'fwh2': 172,
    'fwh3': 183,
    'fwh4': 192,
}


def fwh_surface_radius(grid_file, surface):
    """Return the physical radius of the FWH surface `surface` (fwh1..fwh4),
    read from block 1 (0-indexed) of `grid_file` at the mapped radial index."""
    idx = SURFACE_TO_RADIAL_INDEX[surface]
    g = p3d.Grid(grid_file).set_subzone(1, [idx, 0, 0], [idx, -1, 0]).load()
    r = np.mean(np.sqrt(g.xyz[0][0, :, 0, 0] ** 2 +
                        g.xyz[0][0, :, 0, 1] ** 2))
    return float(r)


def load_mike_pressures(probe_name, num_mikes):
    """Read the mike_{probe_name}_##.dat files back into p[nsteps, num_mikes]."""
    columns = []
    for i in range(num_mikes):
        data = np.loadtxt('mike_%s_%02d.dat' % (probe_name, i + 1))
        columns.append(data[:, 1])
    return np.array(columns).T


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
    return p.parse_args()


def main():
    args = parse_args()
    grid_file = '%s.xyz' % args.prefix

    probe_r = fwh_surface_radius(grid_file, args.surface)
    print('%s: radial index %d -> r/D = %.4f' %
          (args.surface, SURFACE_TO_RADIAL_INDEX[args.surface], probe_r))

    compute_sound(args.prefix, args.x0, args.probe_dt, args.distance,
                  args.theta, args.surface, probe_r)

    p = load_mike_pressures(args.surface, args.num_mikes)
    St, SPL, OASPL = windowed_fft(p, num_windows=args.num_windows,
                                  dt=args.probe_dt, mach_number=args.mach)

    np.savetxt('spl_%s.dat' % args.surface,
               np.column_stack([St, SPL]), fmt='%+.18E')
    with open('oaspl_%s.dat' % args.surface, 'a') as f:
        f.write('%+.6E %+.6E\n' % (args.theta, OASPL))
    print('%s @ theta=%g deg, d=%g: OASPL = %.3f dB' %
          (args.surface, args.theta, args.distance, OASPL))


if __name__ == '__main__':
    main()
