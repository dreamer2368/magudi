import plot3dnasa as p3d
import numpy as np
from numpy.linalg import norm

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
