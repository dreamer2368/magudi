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
