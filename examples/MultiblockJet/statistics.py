"""
Compute Favre-averaged streamwise velocity on the centerline (r/D = 0)
and nozzle lipline (r/D = 0.5) from <input_dir>/<prefix>.mean.q and
write a three-column text file plus two PNG plots suitable for
reproducing the "Current" curves in Figures 8.5(a) and 8.6 of
Vishnampet's thesis.
"""
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from magudi_utils import plot3dnasa as p3d
from magudi_utils.RoundJet import getCenterlineSubarray


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


def centerline_favre_u3(prefix):
    z, Q = getCenterlineSubarray('%s.xyz' % prefix, '%s.mean.q' % prefix)
    return z, Q[:, 3] / Q[:, 0]


def lipline_favre_u3(prefix, r=0.5):
    g = p3d.Grid('%s.xyz' % prefix)
    q = p3d.fromfile('%s.mean.q' % prefix)
    fe = extract_const_r(g, q, r=r)
    ring = fe[0][:, :-1, 0, :]
    rhobar     = np.mean(ring[:, :, 0], axis=1)
    rho_u3_bar = np.mean(ring[:, :, 3], axis=1)
    z = g.set_subzone(0, [0, 0, 0], [0, 0, -1]).load().xyz[0][0, 0, :, 2]
    return z, rho_u3_bar / rhobar


def _plot(x, y, title, path):
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.plot(x, y, 'k-', lw=1.5)
    ax.set_xlabel(r'$(x_3 + x_s)/D$')
    ax.set_ylabel(r'$\tilde{u}_3/U_j$')
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main(input_dir='.', output_dir='.', prefix='MultiblockJet', shift=0.,
         mach_number=1.3, gamma=1.4, r_lipline=0.5):
    in_prefix   = os.path.join(input_dir, prefix)
    out_txt     = os.path.join(output_dir, '%s.favre_u3.txt' % prefix)
    out_center  = os.path.join(output_dir, '%s.favre_u3.centerline.png' % prefix)
    out_lipline = os.path.join(output_dir, '%s.favre_u3.lipline.png' % prefix)

    z_c, u3_c = centerline_favre_u3(in_prefix)
    z_l, u3_l = lipline_favre_u3(in_prefix, r=r_lipline)
    assert np.allclose(z_c, z_l)

    T_ratio = 1. / (1. + 0.5 * (gamma - 1.) * mach_number ** 2)
    U_j = mach_number * np.sqrt(T_ratio)

    x        = z_c + shift
    u3_c_hat = u3_c / U_j
    u3_l_hat = u3_l / U_j

    np.savetxt(
        out_txt, np.column_stack([x, u3_c_hat, u3_l_hat]),
        header='(x3+xs)/D    u_tilde_3/U_j (r/D=0)    u_tilde_3/U_j (r/D=0.5)',
        fmt='%+.15E')

    _plot(x, u3_c_hat, r'Centerline ($r/D = 0$)',           out_center)
    _plot(x, u3_l_hat, r'Nozzle lipline ($r/D = 0.5$)',     out_lipline)


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
    args = parser.parse_args()
    main(input_dir=args.input_dir, output_dir=args.output_dir,
         prefix=args.prefix, shift=args.shift)
