import plot3dnasa as p3d
import numpy as np

def extract_yz(f):
    n = f.get_size()
    args = dict()
    if type(f) == p3d.Function:
        args.update(ncomponents=f.ncomponents)
    fe = type(f)(**args).set_size(
        [n[0][2], n[3][0] + n[0][1] + n[1][0] - 2, 1], True)
    si = {3: [slice(None, None, -1), 0, slice(None)],
          0: [0, slice(None), slice(None)],
          1: [slice(None), 0, slice(None)]}
    so = {3: [slice(None), slice(None, n[3][0] - 1), 0],
          0: [slice(None), slice(n[3][0] - 1, n[3][0] + n[0][1] - 1), 0],
          1: [slice(None), slice(n[3][0] + n[0][1] - 1, None), 0]}
    starts = {3: [1, n[3][1]/2, 0], 0: [n[0][0]/2, 0, 0], 1: [1, n[1][1]/2, 0]}
    ends = {3: [-1, n[3][1]/2, -1], 0: [n[0][0]/2, -1, -1],
            1: [-1, n[1][1]/2, -1]}
    for i in [3, 0, 1]:
        f.set_subzone(i, starts[i], ends[i]).load()
        for j in range(f[0].shape[-1]):
            if type(f) == p3d.Grid:
                fe[0][so[i] + [j]] = f[0][si[i] + [2-j]].T
            else:
                fe[0][so[i] + [j]] = f[0][si[i] + [j]].T
    return fe

def extract_xy(f, k=0):
    n = f.get_size()
    args = dict()
    if type(f) == p3d.Function:
        args.update(ncomponents=f.ncomponents)
    fe = type(f)(**args).set_size(
        [[n[i][0], n[i][1], 1] for i in range(5)], True)
    for i in range(5):
        f.set_subzone(i, [0, 0, k], [-1, -1, k]).load()
        fe[i][:,:,0,:] = f[0][:,:,0,:]
    return fe

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

def compute_sound(prefix, x0, dt, d, theta):
    import os
    import fwhsolver as fwh
    g = p3d.Grid('%s.xyz' % prefix)
    n = g.get_size(0)
    ge = extract_const_r(g, g)
    mikes = fwh.get_mikes(8, x0, d, theta)
    probe_files = ['%s.probe_fwh.%s.dat' % (prefix, s)
                   for s in ['E', 'N', 'W', 'S']]
    nsamples = os.stat(probe_files[0]).st_size / \
               (40 * n[0] * ((n[1] - 1) / 4 + 1))
    solver = fwh.FWHSolver(ge, mikes, nsamples, dt, probe_files=probe_files)
    solver.integrate(chunk_size=50)
    for i, mike in enumerate(mikes):
        with open('mike%02d.dat' % (i + 1), 'w') as f:
            np.savetxt(f, np.array([mike.t, mike.p]).T, fmt='%+.18E')

def extract_axisymmetric(g, f, show_progress=True):
    n = g.get_size()
    z = g.set_subzone(0, [0, 0, 0], [0, 0, -1]).load().xyz[0][0,0,:,2]
    r = np.zeros(n[0][1] / 2 + n[1][0])
    r[1:n[0][1]/2] = np.sqrt(
        np.sum(g.set_subzone(0, [n[0][0]/2, n[0][1]/2 + 1, 0],
                             [n[0][0]/2, -2, 0]).load().xyz[0][0,:,0,0:2] ** 2,
               axis=-1))
    r[n[0][1]/2:] = np.sqrt(
        np.sum(g.set_subzone(2, [0, n[2][1]/2, 0],
                             [-1, n[2][1]/2, 0]).load().xyz[0][:,0,0,0:2] ** 2,
               axis=-1))
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
            if type(f) == p3d.Solution:
                g.set_subzone(i, [0, 0, k], [-1, -1, k]).load()
                theta = np.arctan2(g.xyz[0][:,:,0,1], g.xyz[0][:,:,0,0])
                f[0][:,:,0,2] = f[0][:,:,0,1] * np.cos(theta) + \
                                f[0][:,:,0,2] * np.sin(theta)
                f[0][:,:,0,1] = f[0][:,:,0,3]
                f[0][:,:,0,3] = 0.
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

def centerline_rms_fluctuations(prefix, gamma=1.4):
    z = p3d.fromfile('%s.xyz' % prefix, 0,
                     [0, 0, 0], [0, 0, -1]).xyz[0][0,0,:,2]
    q_m = p3d.fromfile('%s.mean.q' % prefix, 0,
                       [0, 0, 0], [0, 0, -1]).toprimitive(gamma).q[0][0,0,:,:]
    probe_file = '%s.probe_centerline.dat' % prefix
    num_samples = os.stat(probe_file).st_size / (40 * z.size)    
    q = np.reshape(np.fromfile(probe_file), [z.size, 5, num_samples],
                   order='F')
    for i in range(3):
        q[:,i+1,:] /= q[:,0,:]
    q[:,4,:] = (gamma - 1.) * (q[:,4,:] - 0.5 * q[:,0,:] *
                               np.sum(q[:,i+1,:] ** 2 for i in range(3)))
    for i in range(q.shape[-1]):
        q[:,:,i] -= q_m
    a = np.empty([q_m.shape[0], 4])
    a[:,0] = z
    a[:,1] = np.mean(q[:,3,:] * q[:,3,:], axis=-1) # streamwise velocity
    a[:,2] = np.mean(q[:,0,:] * q[:,0,:], axis=-1) # density
    a[:,3] = np.mean(q[:,4,:] * q[:,4,:], axis=-1) # pressure
    a[:,1:] = np.sqrt(a[:,1:])
    np.savetxt('%s.centerline_rms_fluctuations.txt' % prefix, a,
               fmt=a.shape[1] * '%+.15E ')
