#!/usr/bin/env python
import numpy as np
import plot3dnasa as p3d

def mapping_function(x, sigma):
    return np.sinh(sigma * x) / np.sinh(sigma)

def grid(size, eta=0.97):
    from scipy.optimize import fsolve
    x_min = -2.0 * np.pi
    x_max =  2.0 * np.pi
    y_min = -1.0
    y_max =  1.0
    z_min = -2.0 / 3.0 * np.pi
    z_max =  2.0 / 3.0 * np.pi
    dy_min = 0.016
    num_uniform = 7
    g = p3d.Grid().set_size(size, True)
    x = np.linspace(x_min, x_max, g.size[0,0] + 1)[:-1]
    z = np.linspace(z_min, z_max, g.size[0,2] + 1)[:-1]

    xi = np.linspace(-1.0, 1.0, g.size[0,1])
    y = np.sin(eta*xi*np.pi/2.0) / np.sin(eta*np.pi/2.0)
    y += (y_max - y_min)/2.0

    dx, dz = (x_max-x_min)/g.size[0,0], (z_max-z_min)/g.size[0,2]
    dy = y[1:] - y[:-1]
    print ('dx : %.5E, dz : %.5E, dy_min : %.5E, dy_max : %.5E' % (dx,dz,np.min(dy),np.max(dy)) )

    for i in range(x.size):
        g.xyz[0][i,:,:,0] = x[i]
    for j in range(y.size):
        g.xyz[0][:,j,:,1] = y[j]
    for k in range(z.size):
        g.xyz[0][:,:,k,2] = z[k]
    return g

def target_state(g, mach_number=0.5, gamma=1.4):
    s = p3d.Solution().copy_from(g).quiescent(gamma)
    s.q[0][:,:,:,1] = mach_number
    return s.fromprimitive()

def initial_condition(g, external_grid='External/RocFlo-CM.00000000.xyz',
                      external_solution='External/RocFlo-CM.00240000.q',
                      chunk_size=20):
    from scipy.signal import resample
    from scipy.interpolate import interp1d
    x, y, z = g.xyz[0][:,0,0,0], g.xyz[0][0,:,0,1], g.xyz[0][0,0,:,2]
    ge = p3d.Grid(external_grid)
    xe = ge.set_subzone(0, [0, 0, 0], [-2, 0, 0]).load().xyz[0][:,0,0,0]
    ye = ge.set_subzone(0, [0, 0, 0], [0, -1, 0]).load().xyz[0][0,:,0,1]
    ze = ge.set_subzone(0, [0, 0, 0], [0, 0, -2]).load().xyz[0][0,0,:,2]
    if x.size == xe.size and z.size == ze.size:
        st = p3d.fromfile(external_solution, 0, [0, 0, 0], [-2, -1, -2])
    else:
        st = p3d.Solution().set_size([x.size, ye.size, z.size], True)
        ends = [int(d) for d in np.linspace(chunk_size - 1, ye.size - 1,
                                            ye.size / chunk_size)]
        subzones = [([0, ends[i-1] + 1, 0], [-2, ends[i], -2])
                    if i > 0 else ([0, 0, 0], [-2, ends[0], -2])
                    for i in range(len(ends))]
        se = p3d.Solution(external_solution)
        for subzone in subzones:
            se.set_subzone(0, *subzone).load()
            q = np.empty_like(se.q[0], dtype='float64')
            q[:,:,:,:] = se.q[0]
            if x.size != xe.size:
                q = resample(q, x.size, axis=0, window='hann')
            if z.size != ze.size:
                q = resample(q, z.size, axis=2, window='hann')
            st.q[0][:,subzone[0][1]:subzone[1][1]+1,:,:] = q
    s = p3d.Solution().copy_from(g)
    for j in range(5):
        for k in range(z.size):
            for i in range(x.size):
                f = interp1d(ye, st.q[0][i,:,k,j], kind='linear',
                             bounds_error=False, fill_value=st.q[0][i,-1,k,j])
                s.q[0][i,0,k,j] = st.q[0][i,0,k,j]
                s.q[0][i,1:,k,j] = f(y[1:])
    return s

if __name__ == '__main__':
    g = grid([300, 200, 300],eta=0.99)
    g.save('ChannelFlow.xyz')
#    target_state(g).save('BoundaryLayer.target.q')
#    initial_condition(g).save('BoundaryLayer.ic.q')
