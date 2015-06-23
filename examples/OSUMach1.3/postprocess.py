import plot3dnasa as p3d
import numpy as np

def extract_yz(g, f):
    n = g.get_size()
    ge = p3d.Grid().set_size([n[0][2], n[4][0] + n[0][1] + n[2][0] - 2, 1],
                             True)
    if type(f) == p3d.Solution:
        fe = p3d.Solution().copy_from(ge)
    else:
        fe = p3d.Function(ncomponents=f.ncomponents).copy_from(ge)
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
        g.set_subzone(i, starts[i], ends[i]).load()
        f.subzone_from(g).load()
        for j in range(2):
            ge.xyz[0][so[i] + [j]] = g.xyz[0][si[i] + [2-j]].T
        for j in range(fe[0].shape[-1]):
            fe[0][so[i] + [j]] = f[0][si[i] + [j]].T
    return ge, fe

def extract_fwh(g, i=161):
    n = g.get_size()
    g_fwh = p3d.Grid()
    g_fwh.set_size([n[0][2], 4 * (n[1][1] - 1) + 1, 1], True)
    for i in range(1, 5):
        g.set_subzone(i, [i, 0, 0], [i, -2, -1]).load()
        for j in range(3):
            g_fwh.xyz[0][:,(i-1)*(n[1][1]-1):i*(n[1][1]-1),0,j] = \
                g.xyz[0][0,:,:,j].T
    g_fwh.xyz[0][:,-1,:,:] = g_fwh.xyz[0][:,0,:,:]
    return g_fwh
    
