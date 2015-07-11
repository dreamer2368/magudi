import plot3dnasa as p3d
import numpy as np
from numpy.linalg import norm
import os
import itertools

class FWHSolver:
    def __init__(self, g, mikes, nsamples, dt, gamma=1.4):
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
        self.get = None
        self._n = g.get_size(0)

    def get_probe(self, sample_index, probe_files, chunk_size):
        if sample_index >= self._curpos and \
           sample_index < self._curpos + chunk_size:
            return self._buffer[:,:,:,sample_index-self._curpos]
        self._curpos = sample_index - sample_index % chunk_size
        m = (self._buffer.shape[0] - 1) / 4 + 1
        n = 40 * m * self._buffer.shape[1] * chunk_size
        for i, filename in enumerate(probe_files):
            with open(filename) as f:
                f.seek(self._curpos * n / chunk_size)
                self._buffer[i*(m-1):(i+1)*(m-1),:,:,:] = np.reshape(
                    np.fromstring(f.read(n), dtype='<f8'),
                    [m, self._buffer.shape[1], 5, chunk_size],
                    order='F')[:-1,:,:,:]
        return self._buffer[:,:,:,sample_index-self._curpos]

    def integrate(self, probe_files=None, chunk_size=20):
        pbar = None
        try:
            from progressbar import ProgressBar, Percentage, Bar, ETA
            print 'Processing FWH data...'
            pbar = ProgressBar(widgets = [Percentage(), ' ', Bar(
                '=', left = '[', right = ']'), ' ', ETA()],
                               maxval = self.nsamples).start()
        except ImportError:
            pass
        if probe_files:
            self._curpos = self.nsamples
            self._buffer = np.empty([self._n[1], self._n[0], 5, chunk_size],
                                    order='F')
        for i in range(self.nsamples):
            if probe_files:
                q = self.get_probe(i, probe_files, chunk_size)
            else:
                q = self.get(i, *self.get_args)
            q[-1,:,:] = q[0,:,:]
            q[:,:,0] = 1. / q[:,:,0]
            q[:,:,4] = (self.gamma - 1.) * (q[:,:,4] - 0.5 * q[:,:,0] * np.sum(
                    q[:,:,i+1] for i in range(3))) - 1. / self.gamma
            for mike in self.mikes:
                mike.add_contribution(i, q)
            if pbar:
                pbar.update(i)
            if i % 100 == 0:
                for j, mike in enumerate(self.mikes):
                    with open('mike%02d.dat' % (j + 1), 'w') as f:
                        np.savetxt(f, np.array([mike.t, mike.p]).T, fmt='%+.18E')
        if pbar:
            pbar.finish()

    def _compute_normals(self, xyz):
        """Computes the areas of quadrilateral elements and the unit
        normal vector at each cell."""
        a = np.cross(xyz[1:,:-1,0,:] - xyz[:-1,:-1,0,:],
                     xyz[:-1,1:,0,:] - xyz[:-1,:-1,0,:])
        b = np.cross(xyz[:-1,1:,0,:] - xyz[1:,1:,0,:],
                     xyz[1:,:-1,0,:] - xyz[1:,1:,0,:])
        ab = norm(a + b, axis=-1)
        cell_areas = 0.5 * (norm(a, axis=-1) + norm(b, axis=-1))
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
        return norm(xyz - self.xyz, axis=-1).min()

    def max_dist(self, xyz):
        """Distance to the farthest grid point."""
        return norm(xyz - self.xyz, axis=-1).max()

    def set_params(self, xyz, cell_areas, unit_normals, dt, nsteps,
                   nsamples, offset):
        n = np.array(xyz.shape[0:2], 'int64') - 1
        self._allocate(n, nsteps)
        dist = np.empty_like(self.normal_projection)
        for i, s in enumerate(self.slices):
            self.disp[:,:,:,i] = self.xyz - xyz[s + [slice(None)]]
            dist[:,:,i] = norm(self.disp[:,:,:,i], axis=-1)
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
            self.Q[:,:,i,-1] = np.sum(q[s + [j+1]] * self.unit_normals[:,:,j]
                                      for j in range(3))
            self.L[:,:,i,-1] = self.normal_projection[:,:,i] * q[s + [4]] + \
                               q[s + [0]] * self.Q[:,:,i,-1] * \
                np.sum(q[s + [j+1]] * self.disp[:,:,j,i] for j in range(3))
        if sample_index < self.coeff.size - 1:
            return self
        dp = (np.sum((self.L[:,:,:,i] + self.Q[:,:,:,i]) * c
                     for i, c in enumerate(self.coeff)) +
              self.L[:,:,:,self.coeff.size//2] * self.dist_inverse) * \
            self.dist_inverse * self.dp_factor
        for i in range(dp.shape[-1]):
            dp[:,:,i] *= self.cell_areas
        return self._update_signal(sample_index, dp)

    def _update_signal(self, sample_index, dp):
        a = np.trunc(self.advanced_offset + sample_index) - \
            self.signal_offset - 1
        self.p[:-1] += np.histogram(a, range(self.p.size),
                                    (0., self.p.size - 2.),
                                    weights=(1. - self.weights) * dp)[0]
        self.p[1:] += np.histogram(a, range(self.p.size),
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

def monopole(sample_index, g, x0, dt, gamma):
    a = 1e-4
    omega = 2. * np.pi / (16. * dt)
    disp = np.rollaxis(g.xyz[0][:,:,0,:] - x0, axis=1, start=0)
    dist = norm(disp, axis=-1)
    t = sample_index * dt - dist
    n = g.get_size(0)
    q = np.empty([n[1], n[0], 5], order='F')
    q[:,:,0] = 1. + a * omega * np.cos(omega * t) / (4. * np.pi * dist)
    q[:,:,1:4] = a * disp / (4. * np.pi)
    for i in range(3):
        q[:,:,i+1] *= q[:,:,0] * (np.sin(omega * t) / dist ** 3 +
                                  omega * np.cos(omega * t) / dist ** 2)
    q[:,:,4] = q[:,:,0] / (gamma - 1.) - 1. / gamma + 0.5 / q[:,:,0] * np.sum(
        q[:,:,i+1] for i in range(3))
    return q

def extract_yz(f):
    n = f.get_size()
    args = dict()
    if type(f) == p3d.Function:
        args.update(ncomponents=f.ncomponents)
    fe = type(f)(**args).set_size(
        [n[0][2], n[4][0] + n[0][1] + n[2][0] - 2, 1], True)
    si = {4: [slice(None, None, -1), 0, slice(None)],
          0: [0, slice(None), slice(None)],
          2: [slice(None), 0, slice(None)]}
    so = {4: [slice(None), slice(None, n[4][0] - 1), 0],
          0: [slice(None), slice(n[4][0] - 1, n[4][0] + n[0][1] - 1), 0],
          2: [slice(None), slice(n[4][0] + n[0][1] - 1, None), 0]}
    starts = {4: [1, n[4][1]/2, 0], 0: [n[0][0]/2, 0, 0], 2: [1, n[2][1]/2, 0]}
    ends = {4: [-1, n[4][1]/2, -1], 0: [n[0][0]/2, -1, -1],
            2: [-1, n[2][1]/2, -1]}
    for i in [4, 0, 2]:
        f.set_subzone(i, starts[i], ends[i]).load()
        for j in range(f[0].shape[-1]):
            if type(f) == p3d.Grid:
                fe[0][so[i] + [j]] = f[0][si[i] + [2-j]].T
            else:
                fe[0][so[i] + [j]] = f[0][si[i] + [j]].T
    return fe

def extract_const_r(g, f=None, r=0.5, stencil_size=5, show_progress=True):
    from scipy.interpolate import BarycentricInterpolator
    z = g.set_subzone(0, [0, 0, 0], [0, 0, -1]).load().xyz[0][0,0,:,2]
    # Find nearest neighbors and setup Barycentric interpolators.
    g.set_subzone(1, [0, 0, 0], [-1, -2, 0]).load()
    i_nearest = np.array([np.argmin(np.abs(g.xyz[0][:,i,0,0] ** 2 +
                                           g.xyz[0][:,i,0,1] ** 2 - r ** 2))
                          for i in range(g.get_size(0)[1])], dtype=int)
    i_offset = i_nearest.min() - stencil_size / 2
    interpolators = [ None ] * i_nearest.size
    for i, j in enumerate(i_nearest):
        interpolators[i] = BarycentricInterpolator(np.sqrt(
        g.xyz[0][j-stencil_size/2:j+stencil_size/2+1,i,0,0] ** 2 +
        g.xyz[0][j-stencil_size/2:j+stencil_size/2+1,i,0,1] ** 2))
    i_nearest -= i_offset
    ge = p3d.Grid().set_size([z.size, 4 * len(interpolators) + 1, 1], True)
    for k in range(z.size):
        ge.xyz[0][k,:,0,2] = z[k]
    m = 2
    if f:
        args = dict()
        if type(f) == p3d.Function:
            args.update(ncomponents=f.ncomponents)
        fe = type(f)(**args).copy_from(ge)
        m += f.set_subzone(0, [0, 0, 0], [1, 0, 0]).load()[0].shape[-1]
    if show_progress:
        from progressbar import ProgressBar, Percentage, Bar, ETA
        print 'Interpolating on constant-radius surface:'
        p = ProgressBar(widgets = [Percentage(), ' ',
                                   Bar('=', left = '[', right = ']'), ' ',
                                   ETA()], maxval=4 * m).start()
    for i in range(1, 5):
        g.set_subzone(
            i, [i_offset + i_nearest.min() - stencil_size / 2, 0, 0],
            [i_offset + i_nearest.max() + stencil_size / 2, -2, 0]).load()
        for l in range(2): # grid coordinates
            for j, interpolator in enumerate(interpolators):
                interpolator.set_yi(
                    g.xyz[0][i_nearest[j]-stencil_size/2:
                             i_nearest[j]+stencil_size/2+1,j,0,l])
                ge.xyz[0][:,(i-1)*len(interpolators)+j,0,l] = interpolator(r)
            if show_progress:
                p.update((i - 1) * m + l + 1)
        if f:
            f.set_subzone(
                i, [i_offset + i_nearest.min() - stencil_size / 2, 0, 0],
                [i_offset + i_nearest.max() + stencil_size / 2, -2, -1]).load()
        for l in range(m - 2): # solution/function components
            for k in range(z.size):
                for j, interpolator in enumerate(interpolators):
                    interpolator.set_yi(
                        f[0][i_nearest[j]-stencil_size/2:
                             i_nearest[j]+stencil_size/2+1,j,k,l])
                    fe[0][k,(i-1)*len(interpolators)+j,0,l] = interpolator(r)
            if show_progress:
                p.update((i - 1) * m + l + 3)
    if show_progress:
        p.finish()
    ge.xyz[0][:,-1,:,:] = ge.xyz[0][:,0,:,:]
    if f:
        fe[0][:,-1,:,:] = fe[0][:,0,:,:]
        return ge, fe
    return ge

def extract_fwh(g, i_fwh=161):
    n = g.get_size()
    g_fwh = p3d.Grid()
    g_fwh.set_size([n[0][2], 4 * (n[1][1] - 1) + 1, 1], True)
    for i in range(1, 5):
        g.set_subzone(i, [i_fwh, 0, 0], [i_fwh, -2, -1]).load()
        for j in range(3):
            g_fwh.xyz[0][:,(i-1)*(n[1][1]-1):i*(n[1][1]-1),0,j] = \
                g.xyz[0][0,:,:,j].T
    g_fwh.xyz[0][:,-1,:,:] = g_fwh.xyz[0][:,0,:,:]
    return g_fwh

def farfield_sound(g, probe_files=['OSUMach1.3.probe_fwh.%s.dat' % s
                                   for s in ['E', 'N', 'W', 'S']],
                   x0=[0., 0., 1.1], d=94., theta=30., dt=0.048, gamma=1.4):
    g = extract_fwh(g)
    mikes = get_mikes(8, x0, d, theta)
    n = g.get_size(0)
    nsamples = os.stat(probe_files[0]).st_size / \
               (40 * n[0] * ((n[1] - 1) / 4 + 1))
    solver = FWHSolver(g, mikes, nsamples, dt, gamma=gamma)
    # solver.get = monopole
    # solver.get_args = (g, x0, dt, gamma)
    solver.integrate(probe_files)
    return mikes

def windowed_fft(p, num_windows=5, dt=0.048, mach_number=1.3, gamma=1.4):
    import numpy.fft
    n = p.shape[0]
    m = 2 * (n / (num_windows + 1))
    windows = [((int(0.5 * i * m), int(0.5 * i * m) + m))
               for i in range(num_windows)]
    temperature_ratio = 1. / (1. + 0.5 * (gamma - 1.) * mach_number ** 2)
    u_j = mach_number * np.sqrt(temperature_ratio)
    St = numpy.fft.fftfreq(m, d=dt)[1:m/2] / u_j
    y = np.empty([m / 2, num_windows, p.shape[1]])
    window_func = np.blackman(m)
    for j in range(p.shape[1]):
        for i, w in enumerate(windows):
            y[:,i,j] = np.absolute(numpy.fft.fft(
                p[w[0]:w[1],j] * window_func))[:m/2] / (m * window_func.mean())
            y[1:,i,j] *= np.sqrt(2.)
    p_ref = 20.e-6 / 101325. / gamma
    OASPL = 10. * np.log10(np.mean(np.sum(np.mean(
        y[1:] ** 2, axis=1), axis=0)) / p_ref ** 2)
    SPL = 10. * np.log10(np.mean(np.mean(y[1:] ** 2, axis=1), axis=1) /
                         p_ref ** 2)
    SPL += 10. * np.log10(1. / dt / n)
    return St, SPL, OASPL

def extract_axisymmetric(g, f, show_progress=True):
    n = g.get_size()
    z = g.set_subzone(0, [0, 0, 0], [0, 0, -1]).load().xyz[0][0,0,:,2]
    r = np.zeros(n[0][1] / 2 + n[1][0])
    r[1:n[0][1]/2] = g.set_subzone(0, [n[0][0]/2, n[0][1]/2 + 1, 0],
                                  [n[0][0]/2, -2, 0]).load().xyz[0][0,:,0,1]
    r[n[0][1]/2:] = g.set_subzone(2, [0, n[2][1]/2, 0],
                                  [-1, n[2][1]/2, 0]).load().xyz[0][:,0,0,1]
    bins = [r[0] - 0.5 * (r[1] - r[0])] + list(0.5 * (r[:-1] + r[1:])) + \
           [r[-1] + 0.5 * (r[-1] - r[-2])]
    args = dict()
    if type(f) == p3d.Function:
        args.update(ncomponents=f.ncomponents)
    fe = type(f)(**args).set_size([z.size, r.size, 1], True)
    ge = p3d.Grid().copy_from(fe)
    for i in range(r.size):
        ge.xyz[0][:,i,0,1] = r[i]
    m = f.set_subzone(0, [0, 0, 0], [1, 0, 0]).load()[0].shape[-1]
    y = np.empty([np.sum(np.prod(n[:,0:2], axis=1)), m])
    x = np.empty(0)
    for i in range(5):
        xy = g.set_subzone(i, [0, 0, 0], [-1, -1, 0]).load().xyz[0][:,:,0,:]
        x = np.append(x, np.sqrt(xy[:,:,0] ** 2 + xy[:,:,1] ** 2).ravel())
    if show_progress:
        from progressbar import ProgressBar, Percentage, Bar, ETA
        print 'Processing along streamwise direction:'
        p = ProgressBar(widgets = [Percentage(), ' ',
                                   Bar('=', left = '[', right = ']'), ' ',
                                   ETA()], maxval=z.size).start()
    for k in range(z.size):
        j = 0
        for i in range(5):
            f.set_subzone(i, [0, 0, k], [-1, -1, k]).load()
            for l in range(m):
                y[j:j+np.prod(n[i,0:2]),l] = f[0][:,:,0,l].ravel()
            j += np.prod(n[i,0:2])
        for i in range(5):
            fe[0][k,:,0,i] = np.histogram(x, bins, weights=y[:,i])[0] / \
                             np.histogram(x, bins)[0]
        ge.xyz[0][k,:,0,0] = z[k]
        if show_progress:
            p.update(k)
    if show_progress:
        p.finish()
    return ge, fe
