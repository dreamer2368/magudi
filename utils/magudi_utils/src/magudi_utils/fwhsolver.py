#!/usr/bin/env python
# -*- coding: utf-8 -*-

from . import plot3dnasa as p3d
import numpy as np

def get_fromfile(offset, size, n, probe_files):
    n0, n1 = int(n[0]), int(n[1])
    q = np.empty([n1, n0, 5, size], order='F')
    m = (n1 - 1) // 4 + 1
    nbytes = 40 * m * n0 * size
    for i, filename in enumerate(probe_files):
        with open(filename, 'rb') as f:
            f.seek(offset * nbytes // size)
            q[i*(m-1):(i+1)*(m-1),:,:,:] = np.reshape(
                np.frombuffer(f.read(nbytes), dtype='<f8'),
                [m, q.shape[1], 5, size], order='F')[:-1,:,:,:]
    return q

class FWHSolver:
    def __init__(self, g, mikes, nsamples, dt, probe_files=None, gamma=1.4):
        d_min = min((mike.min_dist(g.xyz[0]) for mike in mikes))
        d_max = max((mike.max_dist(g.xyz[0]) for mike in mikes))
        self.nsamples = nsamples
        nsteps = nsamples - (int(np.ceil(d_max / dt)) -
                             int(np.floor(d_min / dt)))
        assert nsteps > 0
        self.mikes = mikes
        self.gamma = gamma
        cell_areas, unit_normals = self._compute_normals(g.xyz[0])
        for mike in mikes:
            mike.set_params(np.rollaxis(g.xyz[0][:,:,0,:], axis=1, start=0),
                            cell_areas, unit_normals, dt, nsteps,
                            nsamples, int(np.ceil(d_max / dt)))
        self.get = get_fromfile
        self.get_args = (g.get_size(0), probe_files)

    def integrate(self, chunk_size=20):
        pbar = None
        try:
            from progressbar import ProgressBar, Percentage, Bar, ETA
            print('Processing FWH data...')
            pbar = ProgressBar(widgets = [Percentage(), ' ', Bar(
                '=', left = '[', right = ']'), ' ', ETA()],
                               maxval = self.nsamples).start()
        except ImportError:
            pass
        for i in range(self.nsamples):
            if i % chunk_size == 0:
                actual = min(chunk_size, self.nsamples - i)
                q = self.get(i, actual, *self.get_args)
                q[-1,:,:,:] = q[0,:,:,:]
                q[:,:,0,:] = 1. / q[:,:,0,:]
                q[:,:,4,:] = (self.gamma - 1.) * (
                    q[:,:,4,:] - 0.5 * q[:,:,0,:] *
                    np.sum(q[:,:,1:4,:], axis=2)) - 1. / self.gamma
            for mike in self.mikes:
                mike.add_contribution(i, q[:,:,:,i%chunk_size])
            if pbar:
                pbar.update(i)
            # if i % 100 == 0:
            #     for j, mike in enumerate(self.mikes):
            #         with open('mike%02d.dat' % (j + 1), 'w') as f:
            #             np.savetxt(f, np.array([mike.t, mike.p]).T, fmt='%+.18E')
        if pbar:
            pbar.finish()

    def _compute_normals(self, xyz):
        """Computes the areas of quadrilateral elements and the unit
        normal vector at each cell."""
        a = np.cross(xyz[1:,:-1,0,:] - xyz[:-1,:-1,0,:],
                     xyz[:-1,1:,0,:] - xyz[:-1,:-1,0,:])
        b = np.cross(xyz[:-1,1:,0,:] - xyz[1:,1:,0,:],
                     xyz[1:,:-1,0,:] - xyz[1:,1:,0,:])
        ab = np.sqrt(np.sum((a + b) ** 2, axis=-1))
        cell_areas = 0.5 * (np.sqrt(np.sum(a ** 2, axis=-1)) + 
                            np.sqrt(np.sum(b ** 2, axis=-1)))
        unit_normals = a + b
        for i in range(unit_normals.shape[-1]):
            unit_normals[:,:,i] /= ab
        return np.asfortranarray(cell_areas.T), np.asfortranarray(
            np.rollaxis(unit_normals, axis=1, start=0))


class Mike:
    def __init__(self, xyz=[0., 0., 0.]):
        self.xyz = xyz
        self.slices = [[slice(None, -1), slice(None, -1)],
                       [slice(None, -1), slice(+1, None)],
                       [slice(+1, None), slice(None, -1)],
                       [slice(+1, None), slice(+1, None)]]
        self.coeff = p3d.fdcoeff([-2, -1, 0, 1, 2], order=1)

    def min_dist(self, xyz):
        """Distance to the closest grid point."""
        return np.sqrt(np.sum((xyz - self.xyz) ** 2, axis=-1)).min()

    def max_dist(self, xyz):
        """Distance to the farthest grid point."""
        return np.sqrt(np.sum((xyz - self.xyz) ** 2, axis=-1)).max()

    def set_params(self, xyz, cell_areas, unit_normals, dt, nsteps,
                   nsamples, offset):
        n = np.array(xyz.shape[0:2], 'int64') - 1
        self._allocate(n, nsteps)
        dist = np.empty_like(self.normal_projection)
        for i, s in enumerate(self.slices):
            self.disp[:,:,:,i] = self.xyz - xyz[tuple(s) + (slice(None),)]
            dist[:,:,i] = np.sqrt(np.sum(self.disp[:,:,:,i] ** 2, axis=-1))
        self.dist_inverse = 1. / dist
        self.advanced_offset = dist / dt - self.coeff.size // 2
        for i in range(3):
            self.disp[:,:,i,:] = self.disp[:,:,i,:] * self.dist_inverse
        self.cell_areas = cell_areas
        self.unit_normals = unit_normals
        for i in range(len(self.slices)):
            self.normal_projection[:,:,i] = np.sum(
                self.unit_normals * self.disp[:,:,:,i], axis=-1)
        self.signal_offset = offset
        self.dp_factor = 1. / (4. * np.pi * len(self.slices))
        self.coeff /= dt
        self.weights = self.advanced_offset - np.trunc(self.advanced_offset)
        self.t = (self.signal_offset + 1 + np.arange(nsteps)) * dt
        return self

    def _allocate(self, n, nsteps):
        self.disp = np.empty([n[0], n[1], 3, len(self.slices)], order='F')
        self.normal_projection = np.empty([n[0], n[1], len(self.slices)],
                                          order='F')
        self.Q = np.empty([n[0], n[1], len(self.slices), self.coeff.size],
                          order='F') # monopole strength
        self.L = np.empty_like(self.Q) # dipole strength
        self.p = np.zeros(nsteps)

    def add_contribution(self, sample_index, q):
        self.Q[:,:,:,:-1] = self.Q[:,:,:,1:]
        self.L[:,:,:,:-1] = self.L[:,:,:,1:]
        for i, s in enumerate(self.slices):
            st = tuple(s)
            self.Q[:,:,i,-1] = np.sum(
                q[st + (slice(1, 4),)] * self.unit_normals, axis=-1)
            self.L[:,:,i,-1] = self.normal_projection[:,:,i] * q[st + (4,)] + \
                               q[st + (0,)] * self.Q[:,:,i,-1] * \
                np.sum(q[st + (slice(1, 4),)] * self.disp[:,:,:,i], axis=-1)
        if sample_index < self.coeff.size - 1:
            return self
        dp = (np.sum((self.L + self.Q) * self.coeff, axis=-1) +
              self.L[:,:,:,self.coeff.size//2] * self.dist_inverse) * \
            self.dist_inverse * self.dp_factor
        for i in range(dp.shape[-1]):
            dp[:,:,i] *= self.cell_areas
        return self._update_signal(sample_index, dp)

    def _update_signal(self, sample_index, dp):
        a = np.trunc(self.advanced_offset + sample_index) - \
            self.signal_offset - 1
        self.p[:-1] += np.histogram(a, list(range(self.p.size)),
                                    (0., self.p.size - 2.),
                                    weights=(1. - self.weights) * dp)[0]
        self.p[1:] += np.histogram(a, list(range(self.p.size)),
                                   (0., self.p.size - 2.),
                                   weights=self.weights * dp)[0]
        return self


def get_mikes(n, x0, d, theta):
    """Returns a list of mikes distributed uniformly on a ring parallel
    to the x-y plane whose center is at a distance `d` * cos(`theta`)
    from `x0`. The first mike is located on the x-z plane, and the mikes
    are ordered in counter-clockwise direction when viewed from
    z=-inf."""
    phi = np.linspace(0., 2. * np.pi, n + 1)[:-1]
    return [Mike([x0[0] + d * np.sin(np.pi * theta / 180.) * np.cos(phi[i]),
                  x0[1] + d * np.sin(np.pi * theta / 180.) * np.sin(phi[i]),
                  x0[2] + d * np.cos(np.pi * theta / 180.)]) for i in range(n)]

def monopole_flow(x, y, t, A0, omega, lam=0., gamma=1.4,
                  c_inf=1., rho_inf=1.):
    """Absolute (rho, u, p) field of a stationary time-harmonic acoustic
    monopole in an infinite quiescent medium (Kim 2012 thesis eq. 5.15,
    with optional exponential decay from the same section).

    Source strength q(tau) = A0 * exp(-lam*tau) * sin(omega*tau) at y;
    the observer field is evaluated using the retarded time
    tau_r = t - |x-y|/c_inf.

    Uses the Kim / Vishnampet nondimensionalization: ambient state is
    (rho_inf, 0, rho_inf * c_inf**2 / gamma) (i.e. p_inf = 1/gamma when
    rho_inf = c_inf = 1).

    Parameters
    ----------
    x : array_like, shape (..., 3)
        Field points.
    y : array_like, shape (3,)
        Source location.
    t : array_like
        Time; must broadcast with r = |x-y| (typically a scalar or a 1-D
        time array combined with x expanded on an extra axis).
    A0, omega, lam : float
        Monopole strength, angular frequency, and exponential decay rate.

    Returns
    -------
    dict with keys 'rho' (shape ...), 'u' (shape (..., 3)), 'p' (shape ...).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    t = np.asarray(t, dtype=float)
    disp = x - y
    r = np.linalg.norm(disp, axis=-1)
    e_hat = disp / r[..., None]
    tau_r = t - r / c_inf
    envelope = A0 * np.exp(-lam * tau_r)
    sin_w = np.sin(omega * tau_r)
    cos_w = np.cos(omega * tau_r)
    q_tau  = envelope * sin_w
    qp_tau = envelope * (omega * cos_w - lam * sin_w)
    p_prime = rho_inf * qp_tau / (4. * np.pi * r)
    rho_prime = p_prime / (c_inf ** 2)
    u_r = q_tau / (4. * np.pi * r ** 2) + qp_tau / (4. * np.pi * r * c_inf)
    u_vec = e_hat * u_r[..., None]
    return dict(rho=rho_inf + rho_prime,
                u=u_vec,
                p=rho_inf * c_inf ** 2 / gamma + p_prime)


def dipole_flow(x, y, t, A0, omega, gamma=1.4,
                c_inf=1., rho_inf=1., axis=1):
    """Absolute (rho, u, p) field of a stationary time-harmonic acoustic
    dipole (Kim 2012 thesis eq. 5.17). Aligned in +axis (default axis=1
    is +y, per Kim section 5.3.2 "maximum sound in y direction").

    Dipole source is f_j(tau) = A0 * sin(omega*tau) * delta_{j, axis}.
    Velocity is derived from the dipole velocity potential; contains
    the standard 1/r, 1/r**2, 1/r**3 Cartesian terms.

    Parameters mirror :func:`monopole_flow` (no lam — Kim uses no
    exponential decay for the dipole test).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    t = np.asarray(t, dtype=float)
    disp = x - y
    r = np.linalg.norm(disp, axis=-1)
    e_hat = disp / r[..., None]
    e_d = e_hat[..., axis]
    tau_r = t - r / c_inf
    f  = A0 * np.sin(omega * tau_r)
    fp = A0 * omega * np.cos(omega * tau_r)
    F  = -A0 * np.cos(omega * tau_r) / omega
    p_prime = -e_d / (4. * np.pi) * (fp / (c_inf * r) + f / (r ** 2))
    rho_prime = p_prime / (c_inf ** 2)
    B = f / (c_inf * r) + F / (r ** 2)
    K = (fp / (c_inf ** 2 * r) + 2. * f / (c_inf * r ** 2)
         + 2. * F / (r ** 3))
    basis = np.zeros(3)
    basis[axis] = 1.
    dA_dx = (basis - e_d[..., None] * e_hat) / r[..., None]
    dB_dx = -e_hat * K[..., None]
    u_prime = (dA_dx * B[..., None] + e_d[..., None] * dB_dx) \
              / (4. * np.pi * rho_inf)
    return dict(rho=rho_inf + rho_prime,
                u=u_prime,
                p=rho_inf * c_inf ** 2 / gamma + p_prime)


def _from_primitive(rho, u, p, gamma):
    """Assemble a q-file chunk of shape [n1, n0, 5, size] (Fortran-ordered,
    conservative variables) from (rho, u, p) fields of shapes
    (n1, n0, size) and (n1, n0, size, 3)."""
    n1, n0, size = rho.shape
    q = np.empty([n1, n0, 5, size], order='F')
    q[:, :, 0, :] = rho
    q[:, :, 1:4, :] = np.moveaxis(rho[..., None] * u, -1, -2)
    q[:, :, 4, :] = p / (gamma - 1.) + 0.5 * rho * np.sum(u * u, axis=-1)
    return q


def get_monopole(offset, size, xyz, dt, y, A0, omega, lam=0., gamma=1.4,
                 c_inf=1., rho_inf=1.):
    """FWH-surface chunk adapter for a monopole. Matches the FWHSolver.get
    contract ``get(offset, size, *get_args) -> q[n1, n0, 5, size]``.

    Use as::

        solver.get      = get_monopole
        solver.get_args = (xyz, dt, y, A0, omega, lam)
    """
    xyz = np.asarray(xyz, dtype=float)
    t = dt * np.arange(offset, offset + size, dtype=float)
    flow = monopole_flow(xyz[:, :, None, :], y, t, A0, omega, lam, gamma,
                         c_inf, rho_inf)
    return _from_primitive(flow['rho'], flow['u'], flow['p'], gamma)


def get_dipole(offset, size, xyz, dt, y, A0, omega, gamma=1.4,
               c_inf=1., rho_inf=1., axis=1):
    """FWH-surface chunk adapter for a +axis-aligned dipole (default +y).
    Same output layout as :func:`get_monopole`."""
    xyz = np.asarray(xyz, dtype=float)
    t = dt * np.arange(offset, offset + size, dtype=float)
    flow = dipole_flow(xyz[:, :, None, :], y, t, A0, omega, gamma,
                       c_inf, rho_inf, axis)
    return _from_primitive(flow['rho'], flow['u'], flow['p'], gamma)


def make_cylindrical_grid(n_axial, n_azimuthal, radius,
                          z_range=(-9., 34.)):
    """Constant-radius cylindrical FWH surface as a single-block plot3dnasa
    Grid of shape (n_axial, n_azimuthal + 1, 1). The azimuthal direction
    wraps (last row = first row), matching the extract_const_r convention
    in examples/MultiblockJet/postprocess.py. Cylinder axis is +z."""
    z = np.linspace(z_range[0], z_range[1], n_axial)
    theta = np.linspace(0., 2. * np.pi, n_azimuthal + 1)
    g = p3d.Grid().set_size([n_axial, n_azimuthal + 1, 1], True)
    xyz = g.xyz[0]
    xyz[:, :, 0, 0] = radius * np.cos(theta)[None, :]
    xyz[:, :, 0, 1] = radius * np.sin(theta)[None, :]
    xyz[:, :, 0, 2] = z[:, None]
    return g


def run_fwh_reference(kind, source_args, mikes, radius, n_axial, n_azimuthal,
                      dt, nsamples, z_range=(-9., 34.), chunk_size=50,
                      gamma=1.4):
    """Build a synthetic cylindrical FWH surface, wire the analytic
    monopole/dipole generator into FWHSolver, and populate the mikes.

    kind         : 'monopole' or 'dipole'.
    source_args  : positional args after (offset, size, xyz, dt) for the
                   chosen adapter. Monopole: (y, A0, omega[, lam]).
                   Dipole:   (y, A0, omega[, gamma[, c_inf, rho_inf, axis]]).
    Returns the input `mikes` list; each mike has `.t` and `.p` populated.
    """
    g = make_cylindrical_grid(n_axial, n_azimuthal, radius, z_range=z_range)
    solver = FWHSolver(g, mikes, nsamples, dt, probe_files=None, gamma=gamma)
    xyz_surface = np.rollaxis(g.xyz[0][:, :, 0, :], axis=1, start=0)
    get_fn = {'monopole': get_monopole, 'dipole': get_dipole}[kind]
    solver.get = get_fn
    solver.get_args = (xyz_surface, dt) + tuple(source_args)
    solver.integrate(chunk_size=chunk_size)
    return mikes


def windowed_fft(p, num_windows=5, dt=1.2e-3 * 35, window_type='blackman'):
    import numpy.fft
    from scipy.signal import get_window
    n = p.shape[0]
    m = 2 * (n // (num_windows + 1))
    windows = [((int(0.5 * i * m), int(0.5 * i * m) + m))
               for i in range(num_windows)]
    y = np.empty([(m + 1) // 2, num_windows, p.shape[1]])
    if window_type:
        window_func = get_window(window_type, m)
    else:
        window_func = np.ones(m)
    for j in range(p.shape[1]):
        for i, w in enumerate(windows):
            y[:,i,j] = np.absolute(numpy.fft.fft(
                p[w[0]:w[1],j] * window_func))[:(m+1)//2] / window_func.sum()
            y[1:,i,j] *= np.sqrt(2.)
    p_hat = np.sqrt(np.mean(np.mean(y ** 2, axis=1), axis=1))
    return numpy.fft.fftfreq(m, d=dt)[:p_hat.size], p_hat
