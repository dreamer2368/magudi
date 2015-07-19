import plot3dnasa as p3d
import numpy as np
from numpy.linalg import norm
import os
import itertools

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
