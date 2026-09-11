"""Tests for magudi_utils.fwhsolver.

Sanity tests for the reference-flow generators, plus round-trip FWH tests
reproducing Kim (2012) thesis figures 5.6-5.9.

Geometry convention:
- Kim uses cylinder axis = +x and defines the polar angle phi w.r.t. +x.
- Our FWH surface (matching MultiblockJet) uses cylinder axis = +z; the
  dipole is aligned in +y. To reproduce Kim's polar sweep we place mikes
  in the y-z plane at various polar angles phi measured from +z (our
  cylinder axis), so phi=0 is on-axis (+z), phi=pi/2 is +y (perpendicular
  to axis, aligned with dipole peak), phi=pi is -z, phi=3*pi/2 is -y.

Kim's source parameters (section 5.3):
    A0    = 1e-4
    omega = 2*pi / (dt*PPW), dt=0.18, PPW=16   -> omega ~= 2.18, f ~= 0.35
    d     = 40 (observer distance)
    lam   = 0.01 (decaying monopole only)

FWH surface: constant-radius cylinder along +z with axial extent [-9, 34],
matching examples/MultiblockJet/compute_sound.py. The four MultiblockJet
radii (block-1 radial indices 155/172/183/192 for grid(60, 196, 132) in
examples/MultiblockJet/config.py) are hard-coded below.

Set MAGUDI_TEST_PLOTS=1 to also write PNGs for visual comparison with Kim.
"""
import os
import numpy as np
import pytest

from magudi_utils import fwhsolver as fwh
from magudi_utils import plot3dnasa as p3d


GAMMA = 1.4
P_INF = 1.0 / GAMMA

KIM_A0    = 1e-4
KIM_PPW   = 16
KIM_DT    = 0.18
KIM_OMEGA = 2. * np.pi / (KIM_DT * KIM_PPW)          # ~ 2.1817
KIM_LAM   = 0.01
KIM_D     = 40.0
Y_SRC     = np.zeros(3)

FWH_RADII = {'fwh1': 1.513909, 'fwh2': 2.045137,
             'fwh3': 2.527544, 'fwh4': 3.031904}
DEFAULT_RADIUS = FWH_RADII['fwh2']
Z_RANGE       = (-9., 34.)
N_AXIAL_TEST  = 512
N_AZIM_TEST   = 128

FIG_DIR = os.path.join(os.path.dirname(__file__), 'figures')


def _plots_enabled():
    return bool(os.environ.get('MAGUDI_TEST_PLOTS'))


def _ensure_fig_dir():
    if _plots_enabled():
        os.makedirs(FIG_DIR, exist_ok=True)


def _mike_on_yz(phi, d=KIM_D):
    """Mike at polar angle `phi` from +z (our cyl axis) in the y-z plane."""
    return fwh.Mike([0., d * np.sin(phi), d * np.cos(phi)])


# ------------------------------------------------------------------
# Sanity tests.
# ------------------------------------------------------------------

def test_monopole_flow_shapes_and_ambient():
    """With A0=0 the analytic field is exactly the quiescent ambient state,
    and the return shapes follow the documented broadcasting contract."""
    x = np.array([0., 0., 10.])
    y = np.zeros(3)
    flow = fwh.monopole_flow(x, y, t=0., A0=0., omega=1.0)
    assert flow['rho'].shape == ()
    assert flow['u'].shape == (3,)
    assert flow['p'].shape == ()
    assert flow['rho'] == pytest.approx(1.0)
    np.testing.assert_allclose(flow['u'], 0.0)
    assert flow['p'] == pytest.approx(P_INF)


def test_dipole_flow_directivity_signs():
    """A +y-aligned dipole radiates zero p' in the equatorial plane and
    opposite-sign p' at the +y and -y poles at the same instant."""
    y = np.zeros(3)
    r = 10.
    t = 0.7
    kw = dict(A0=1e-4, omega=2.0)
    p_plus  = fwh.dipole_flow(np.array([0.,  r, 0.]), y, t, **kw)['p'] - P_INF
    p_minus = fwh.dipole_flow(np.array([0., -r, 0.]), y, t, **kw)['p'] - P_INF
    p_eq    = fwh.dipole_flow(np.array([r,  0., 0.]), y, t, **kw)['p'] - P_INF
    assert p_plus * p_minus < 0.
    assert abs(p_eq) < 1e-14


def test_get_monopole_shape_and_ambient_offset():
    """The FWH-surface adapter returns a Fortran-ordered [n1, n0, 5, size]
    chunk and, with A0=0, holds the quiescent ambient state."""
    n1, n0, size = 5, 7, 3
    theta = np.linspace(0., 2. * np.pi, n1, endpoint=False)
    xax = np.linspace(-2., 2., n0)
    R = 3.
    xyz = np.stack(np.broadcast_arrays(
        xax[None, :],
        R * np.cos(theta)[:, None],
        R * np.sin(theta)[:, None],
    ), axis=-1)
    q = fwh.get_monopole(offset=0, size=size, xyz=xyz, dt=0.1,
                         y=np.zeros(3), A0=0., omega=1.0)
    assert q.shape == (n1, n0, 5, size)
    assert q.flags['F_CONTIGUOUS']
    np.testing.assert_allclose(q[:, :, 0, :], 1.0)
    np.testing.assert_allclose(q[:, :, 1:4, :], 0.0, atol=1e-15)
    np.testing.assert_allclose(q[:, :, 4, :], P_INF / (GAMMA - 1.))


# ------------------------------------------------------------------
# Geometry: verify _compute_normals returns outward-facing normals.
# ------------------------------------------------------------------

def test_compute_normal_outward():
    """For a cylindrical FWH surface centered on the +z axis the unit
    normal at each cell must point radially outward, i.e., in the +r
    direction (cos theta_center, sin theta_center, 0). We compare the
    normal from FWHSolver._compute_normals with the geometric outward
    direction at the cell center and require unit_normal . outward > 0
    (a positive dot product; magnitude ~ 1 for a well-formed cylinder)."""
    n_axial, n_azim = 20, 32
    radius = 2.0
    g = fwh.make_cylindrical_grid(n_axial, n_azim, radius,
                                  z_range=(-4., 4.))
    # FWHSolver construction requires at least one mike; put it far away.
    dummy_mike = fwh.Mike([100., 0., 0.])
    solver = fwh.FWHSolver(g, [dummy_mike], nsamples=10, dt=1.,
                           probe_files=None)
    cell_areas, unit_normals = solver._compute_normals(g.xyz[0])
    # unit_normals shape after np.rollaxis: [n_azim, n_axial-1, 3]
    xyz = g.xyz[0]  # [n_axial, n_azim+1, 1, 3]
    # Cell (i_azim, i_axial) is bounded by the 4 corner grid points; take
    # the geometric mean of the 4 (x, y) to get the cell-center azimuthal.
    corners = np.stack([
        xyz[:-1, :-1, 0, :2],
        xyz[+1:, :-1, 0, :2],
        xyz[:-1, +1:, 0, :2],
        xyz[+1:, +1:, 0, :2],
    ]).mean(axis=0)  # [n_axial-1, n_azim, 2]
    r_center = np.linalg.norm(corners, axis=-1)  # [n_axial-1, n_azim]
    outward_xy = corners / r_center[..., None]   # unit outward radial
    outward = np.concatenate([outward_xy,
                              np.zeros(outward_xy.shape[:-1] + (1,))],
                             axis=-1)  # [n_axial-1, n_azim, 3]
    outward = np.rollaxis(outward, axis=1, start=0)  # [n_azim, n_axial-1, 3]
    dot = np.sum(unit_normals * outward, axis=-1)
    # All dot products must be positive (outward-facing) and close to +1
    # for a well-resolved cylinder.
    assert np.all(dot > 0.), (
        '_compute_normals returned inward-facing normals; min dot = %.4f'
        % float(dot.min()))
    np.testing.assert_allclose(dot, 1.0, atol=1e-3)


@pytest.fixture(scope='module')
def multiblockjet_fwh_grid(tmp_path_factory):
    """Factory returning the MultiblockJet FWH surface at a given radius.

    Builds the production MultiblockJet grid once, writes it to a temp
    PLOT3D file (extract_const_r requires a file-backed Grid, since it
    calls .load()), and returns a callable `extract(radius)` that runs
    the actual examples/MultiblockJet/postprocess.extract_const_r on it.
    Skipped if the examples/MultiblockJet sources are not importable."""
    import sys
    mbj_dir = os.path.abspath(os.path.join(
        os.path.dirname(__file__), '..', '..', '..',
        'examples', 'MultiblockJet'))
    if not os.path.isdir(mbj_dir):
        pytest.skip('examples/MultiblockJet not found at %s' % mbj_dir)
    if mbj_dir not in sys.path:
        sys.path.insert(0, mbj_dir)
    try:
        from config import grid
        from postprocess import extract_const_r
    except ImportError as e:
        pytest.skip('cannot import MultiblockJet config/postprocess: %s' % e)
    g_full = grid(60, 196, 132, a_inner=0.24, p_inner=1.08634735266)
    grid_path = str(tmp_path_factory.mktemp('mbj') / 'MultiblockJet.xyz')
    g_full.save(grid_path)
    def extract(radius):
        g = p3d.Grid(grid_path)
        return extract_const_r(g, g, r=radius)
    return extract


def test_compute_normal_outward_multiblockjet(multiblockjet_fwh_grid):
    """Same outward-normal check as test_compute_normal_outward, but on the
    actual FWH surface produced by examples/MultiblockJet/postprocess.
    extract_const_r (stitched from blocks 1..4 of the production grid)."""
    ge = multiblockjet_fwh_grid(FWH_RADII['fwh2'])
    dummy_mike = fwh.Mike([100., 0., 0.])
    solver = fwh.FWHSolver(ge, [dummy_mike], nsamples=200, dt=1.0,
                           probe_files=None)
    cell_areas, unit_normals = solver._compute_normals(ge.xyz[0])
    xyz = ge.xyz[0]  # [n_axial, n_azim_full, 1, 3]
    corners = np.stack([
        xyz[:-1, :-1, 0, :2],
        xyz[+1:, :-1, 0, :2],
        xyz[:-1, +1:, 0, :2],
        xyz[+1:, +1:, 0, :2],
    ]).mean(axis=0)  # [n_axial-1, n_azim_cells, 2]
    r_center = np.linalg.norm(corners, axis=-1)
    outward_xy = corners / r_center[..., None]
    outward = np.concatenate([outward_xy,
                              np.zeros(outward_xy.shape[:-1] + (1,))],
                             axis=-1)
    outward = np.rollaxis(outward, axis=1, start=0)
    dot = np.sum(unit_normals * outward, axis=-1)
    assert np.all(dot > 0.), (
        'MultiblockJet FWH grid normals face inward; min dot = %.4f'
        % float(dot.min()))
    np.testing.assert_allclose(dot, 1.0, atol=1e-2)


# ------------------------------------------------------------------
# Round-trip tests (Kim thesis figures 5.6-5.9).
# ------------------------------------------------------------------

def _analytic_p_at(mike, kind, t, **src_kwargs):
    x = np.asarray(mike.xyz)
    flow_fn = {'monopole': fwh.monopole_flow,
               'dipole':   fwh.dipole_flow}[kind]
    return flow_fn(x, Y_SRC, t, **src_kwargs)['p'] - P_INF


def _rel_linf_error(observed, reference):
    return float(np.max(np.abs(observed - reference)) /
                 np.max(np.abs(reference)))


def test_fwh_monopole_time_history_phi30():
    """Kim fig 5.6(a): FWH-reconstructed p'(t) at polar angle phi=30 deg
    from cyl axis should match the analytic monopole at that observer."""
    mike = _mike_on_yz(np.pi / 6.)
    nsamples = 900
    src = (Y_SRC, KIM_A0, KIM_OMEGA)
    fwh.run_fwh_reference('monopole', src, [mike],
                          radius=DEFAULT_RADIUS,
                          n_axial=N_AXIAL_TEST, n_azimuthal=N_AZIM_TEST,
                          dt=KIM_DT, nsamples=nsamples)
    p_ref = _analytic_p_at(mike, 'monopole', mike.t,
                           A0=KIM_A0, omega=KIM_OMEGA)
    if _plots_enabled():
        _ensure_fig_dir()
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.plot(mike.t, mike.p, 'r-', lw=1.2, label='Current (FWH)')
        ax.plot(mike.t, p_ref, 'k--', lw=1.0, label='Monopole (analytic)')
        ax.set_xlim(100., 115.)
        ax.set_xlabel(r'$t$')
        ax.set_ylabel(r"$p'(t)$")
        ax.legend(fontsize=8, frameon=False)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG_DIR, 'kim_5_6a.png'), dpi=150)
        plt.close(fig)
    # Compare only on Kim's plot window; earlier samples include the
    # source-startup transient that Kim also skips.
    m = (mike.t >= 100.) & (mike.t <= 115.)
    err = _rel_linf_error(mike.p[m], p_ref[m])
    assert err < 0.05, 'L-inf relative error %.4f exceeds 5%%' % err


def test_fwh_monopole_spectrum_phi30():
    """Kim fig 5.6(b): FFT of the FWH-reconstructed signal should peak at
    the source frequency f = omega/(2*pi)."""
    mike = _mike_on_yz(np.pi / 6.)
    nsamples = 900
    src = (Y_SRC, KIM_A0, KIM_OMEGA)
    fwh.run_fwh_reference('monopole', src, [mike],
                          radius=DEFAULT_RADIUS,
                          n_axial=N_AXIAL_TEST, n_azimuthal=N_AZIM_TEST,
                          dt=KIM_DT, nsamples=nsamples)
    p = mike.p - np.mean(mike.p)
    n = p.size
    window = np.blackman(n)
    P_hat = np.abs(np.fft.rfft(p * window)) / window.sum()
    P_hat[1:] *= 2.0
    freqs = np.fft.rfftfreq(n, d=KIM_DT)
    peak_idx = int(np.argmax(P_hat))
    f_expected = KIM_OMEGA / (2. * np.pi)
    df = freqs[1] - freqs[0]
    assert abs(freqs[peak_idx] - f_expected) < 1.5 * df
    if _plots_enabled():
        _ensure_fig_dir()
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.plot(freqs, P_hat, 'r-', lw=1.2)
        ax.set_xlim(0., 2.7)
        ax.set_xlabel(r'$f$')
        ax.set_ylabel(r'$|\hat{p}|$')
        fig.tight_layout()
        fig.savefig(os.path.join(FIG_DIR, 'kim_5_6b.png'), dpi=150)
        plt.close(fig)


def test_fwh_decaying_monopole_time_history_phi30():
    """Kim fig 5.8: decaying monopole (lam=0.01) reconstructed at phi=30 deg."""
    mike = _mike_on_yz(np.pi / 6.)
    nsamples = 1000
    src = (Y_SRC, KIM_A0, KIM_OMEGA, KIM_LAM)
    fwh.run_fwh_reference('monopole', src, [mike],
                          radius=DEFAULT_RADIUS,
                          n_axial=N_AXIAL_TEST, n_azimuthal=N_AZIM_TEST,
                          dt=KIM_DT, nsamples=nsamples)
    p_ref = _analytic_p_at(mike, 'monopole', mike.t,
                           A0=KIM_A0, omega=KIM_OMEGA, lam=KIM_LAM)
    if _plots_enabled():
        _ensure_fig_dir()
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.plot(mike.t, mike.p, 'r-', lw=1.2, label='Current (FWH)')
        ax.plot(mike.t, p_ref, 'k--', lw=1.0,
                label='Decaying monopole (analytic)')
        ax.set_xlim(70., 120.)
        ax.set_xlabel(r'$t$')
        ax.set_ylabel(r"$p'(t)$")
        ax.legend(fontsize=8, frameon=False)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG_DIR, 'kim_5_8.png'), dpi=150)
        plt.close(fig)
    m = (mike.t >= 70.) & (mike.t <= 120.)
    err = _rel_linf_error(mike.p[m], p_ref[m])
    assert err < 0.05, 'L-inf relative error %.4f exceeds 5%%' % err


def test_fwh_dipole_time_history_phi60():
    """Companion to test_fwh_monopole_time_history_phi30 for the +y-aligned
    dipole. Mike at polar angle phi=60 deg from cyl axis in the y-z plane;
    same Kim source parameters and comparison window."""
    mike = _mike_on_yz(np.pi / 3.)
    nsamples = 900
    src = (Y_SRC, KIM_A0, KIM_OMEGA)
    fwh.run_fwh_reference('dipole', src, [mike],
                          radius=DEFAULT_RADIUS,
                          n_axial=N_AXIAL_TEST, n_azimuthal=N_AZIM_TEST,
                          dt=KIM_DT, nsamples=nsamples)
    p_ref = _analytic_p_at(mike, 'dipole', mike.t,
                           A0=KIM_A0, omega=KIM_OMEGA)
    if _plots_enabled():
        _ensure_fig_dir()
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.plot(mike.t, mike.p, 'r-', lw=1.2, label='Current (FWH)')
        ax.plot(mike.t, p_ref, 'k--', lw=1.0, label='Dipole (analytic)')
        ax.set_xlim(100., 115.)
        ax.set_xlabel(r'$t$')
        ax.set_ylabel(r"$p'(t)$")
        ax.legend(fontsize=8, frameon=False)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG_DIR, 'dipole_phi60.png'), dpi=150)
        plt.close(fig)
    m = (mike.t >= 100.) & (mike.t <= 115.)
    err = _rel_linf_error(mike.p[m], p_ref[m])
    assert err < 0.05, 'L-inf relative error %.4f exceeds 5%%' % err


def _ring_mikes_polar(n=24, d=KIM_D):
    """Ring of n mikes at various polar angles phi from +z (cyl axis),
    uniformly distributed in the y-z plane at distance d from origin.
    Reproduces Kim's polar-angle sweep (fig 5.7 / 5.9)."""
    phi = np.linspace(0., 2. * np.pi, n, endpoint=False)
    return [_mike_on_yz(p, d) for p in phi], phi


def _mike_rms(mike):
    p = mike.p - np.mean(mike.p)
    return float(np.sqrt(np.mean(p ** 2)))


def _analytic_rms(mike, kind, **src_kwargs):
    p = _analytic_p_at(mike, kind, mike.t, **src_kwargs)
    p = p - np.mean(p)
    return float(np.sqrt(np.mean(p ** 2)))


def test_fwh_monopole_directivity():
    """Kim fig 5.7: monopole directivity is uniform. FWH open-cylinder
    errors are confined to a narrow cone around the cyl axis (phi ~ 0, pi);
    assert uniformity over the equatorial band (30 deg <= phi <= 150 deg)."""
    n_mikes = 24
    mikes, phi = _ring_mikes_polar(n_mikes)
    nsamples = 700
    src = (Y_SRC, KIM_A0, KIM_OMEGA)
    fwh.run_fwh_reference('monopole', src, mikes,
                          radius=DEFAULT_RADIUS,
                          n_axial=N_AXIAL_TEST, n_azimuthal=N_AZIM_TEST,
                          dt=KIM_DT, nsamples=nsamples)
    rms = np.array([_mike_rms(m) for m in mikes])
    # rms_analytic = np.array([_analytic_rms(m, 'monopole',
    #                                        A0=KIM_A0, omega=KIM_OMEGA)
    #                          for m in mikes])
    predicted = KIM_A0 * KIM_OMEGA / (4. * np.pi * KIM_D) / np.sqrt(2.)
    if _plots_enabled():
        _ensure_fig_dir()
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(4.5, 4.5))
        ax = fig.add_subplot(111, projection='polar')
        ax.plot(phi, rms ** 2, 'r-s', lw=1.0, ms=4, label='Current (FWH)')
        # ax.plot(phi, rms_analytic ** 2, 'b--', lw=1.0,
        #         label='Analytic p(mike, t) RMS')
        ax.plot(phi, np.full_like(phi, predicted ** 2),
                'k--o', lw=0.8, ms=3, label='Monopole (analytic)')
        ax.set_title('Kim fig 5.7 (monopole)')
        ax.legend(fontsize=7, loc='lower left', bbox_to_anchor=(1.05, 0.))
        fig.tight_layout()
        fig.savefig(os.path.join(FIG_DIR, 'kim_5_7.png'), dpi=150)
        plt.close(fig)
    # Equatorial band (avoiding the cyl-axis cone where open-surface errors
    # dominate, per Kim p. 70): 30 deg <= phi <= 150 deg (& 210..330).
    eq = ((phi >= np.pi / 6.) & (phi <= 5 * np.pi / 6.)) | \
         ((phi >= 7 * np.pi / 6.) & (phi <= 11 * np.pi / 6.))
    rms_eq = rms[eq]
    assert rms_eq.max() / rms_eq.min() - 1. < 0.10
    assert abs(rms_eq.mean() / predicted - 1.) < 0.05


def test_fwh_dipole_directivity():
    """Kim fig 5.9: dipole (aligned +y) directivity is |sin(phi)| in polar
    angle from cyl axis (+z) — figure-8 with lobes at phi = pi/2 (+y) and
    phi = 3*pi/2 (-y), nulls at phi = 0, pi (on cyl axis)."""
    n_mikes = 24
    mikes, phi = _ring_mikes_polar(n_mikes)
    nsamples = 700
    src = (Y_SRC, KIM_A0, KIM_OMEGA)
    fwh.run_fwh_reference('dipole', src, mikes,
                          radius=DEFAULT_RADIUS,
                          n_axial=N_AXIAL_TEST, n_azimuthal=N_AZIM_TEST,
                          dt=KIM_DT, nsamples=nsamples)
    rms = np.array([_mike_rms(m) for m in mikes])
    # rms_analytic = np.array([_analytic_rms(m, 'dipole',
    #                                        A0=KIM_A0, omega=KIM_OMEGA)
    #                          for m in mikes])
    peak_amp = KIM_A0 * KIM_OMEGA / (4. * np.pi * KIM_D) / np.sqrt(2.)
    predicted = peak_amp * np.abs(np.sin(phi))
    if _plots_enabled():
        _ensure_fig_dir()
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(4.5, 4.5))
        ax = fig.add_subplot(111, projection='polar')
        ax.plot(phi, rms ** 2, 'r-s', lw=1.0, ms=4, label='Current (FWH)')
        # ax.plot(phi, rms_analytic ** 2, 'b--', lw=1.0,
        #         label='Analytic p(mike, t) RMS')
        ax.plot(phi, predicted ** 2, 'k--o', lw=0.8, ms=3,
                label='Dipole (analytic)')
        ax.set_title('Kim fig 5.9 (dipole, +y-aligned)')
        ax.legend(fontsize=7, loc='lower left', bbox_to_anchor=(1.05, 0.))
        fig.tight_layout()
        fig.savefig(os.path.join(FIG_DIR, 'kim_5_9.png'), dpi=150)
        plt.close(fig)
    peak_meas = rms.max()
    assert abs(peak_meas / peak_amp - 1.) < 0.10
    # Nulls (phi = 0, pi) are much smaller than peak:
    for i in (0, n_mikes // 2):
        assert rms[i] < 0.10 * peak_meas
    np.testing.assert_allclose(rms, predicted, atol=0.10 * peak_amp)
