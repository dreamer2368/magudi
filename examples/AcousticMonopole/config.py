#!/usr/bin/env python
import numpy as np
import plot3dnasa as p3d

def grid(size):
    x_min = -14.
    x_max =  14.
    y_min = -14.
    y_max =  14.
    g = p3d.Grid().set_size(size, True)
    x = np.linspace(x_min, x_max, g.size[0,0])
    y = np.linspace(y_min, y_max, g.size[0,1])
    g.xyz[0][:,:,0,:2] = np.transpose(np.meshgrid(x, y))
    return g

def target_mollifier(g):
    x_min =  -1.
    x_max =   1.
    y_min = -10.
    y_max =  10.
    f = p3d.Function().copy_from(g)
    f.f[0].fill(1.)
    n = f.get_size(0)
    for i in range(n[0]):
        f.f[0][i,:,0,0] *= p3d.tanh_support(
            g.xyz[0][i,:,0,1], y_min, y_max, 40., 0.2)
    for j in range(n[1]):
        f.f[0][:,j,0,0] *= p3d.cubic_bspline_support(
            g.xyz[0][:,j,0,0], x_min, x_max)
    imin, imax = p3d.find_extents(g.xyz[0][:,0,0,0], x_min, x_max)
    jmin, jmax = p3d.find_extents(g.xyz[0][0,:,0,1], y_min, y_max)
    print ('  {:<20} {:<21} {:>4d} {:>7d}' + 6 * ' {:>4d}').format(
        'targetRegion', 'COST_TARGET', 1, 0, imin, imax, jmin, jmax, 1, -1)
    return f

def control_mollifier(g):
    x_min =  1.
    x_max =  5.
    y_min = -2.
    y_max =  2.
    f = p3d.Function().copy_from(g)
    f.f[0].fill(1.)
    n = f.get_size(0)
    for j in range(n[1]):
        f.f[0][:,j,0,0] *= p3d.cubic_bspline_support(
            g.xyz[0][:,j,0,0], x_min, x_max)
    for i in range(n[0]):
        f.f[0][i,:,0,0] *= p3d.cubic_bspline_support(
            g.xyz[0][i,:,0,1], y_min, y_max)
    imin, imax = p3d.find_extents(g.xyz[0][:,0,0,0], x_min, x_max)
    jmin, jmax = p3d.find_extents(g.xyz[0][0,:,0,1], y_min, y_max)
    print ('  {:<20} {:<21} {:>4d} {:>7d}' + 6 * ' {:>4d}').format(
        'controlRegion', 'ACTUATOR', 1, 0, imin, imax, jmin, jmax, 1, -1)
    return f

def mean_pressure(s):
    f = p3d.Function().copy_from(s)
    f.f[0][:,:,:,0] = s.toprimitive().q[0][:,:,:,4]
    return f

def random_solution(g,time=0.0,timestep=0):

    gamma = 1.4
    s = p3d.Solution().copy_from(g).quiescent(gamma)
    n = s.get_size(0)
    s.q[0][:,:,0,0] = ( 2.0 * np.random.rand(n[0],n[1]) - 1.0 ) * 1.0e-2 + 1.0
    s.q[0][:,:,0,1] = ( 2.0 * np.random.rand(n[0],n[1]) - 1.0 ) * 1.0e-3 + 0.0
    s.q[0][:,:,0,2] = ( 2.0 * np.random.rand(n[0],n[1]) - 1.0 ) * 1.0e-3 + 0.0
    s.q[0][:,:,0,4] = ( 2.0 * np.random.rand(n[0],n[1]) - 1.0 ) * 1.0e-2 + 1./(gamma-1.)
    s.q[0][:,:,0,4] *= (gamma-1.)/gamma * s.q[0][:,:,0,0]

    s.time = time
    s._format.aux_header[0] = timestep

    return s.fromprimitive(gamma)

if __name__ == '__main__':
    g = grid([201, 201])
    g.save('AcousticMonopole.xyz')

    dt = 5.0e-2
    Nt = 100
    for k in range(2):
        random_solution(g,k*Nt*dt,k*Nt).save('AcousticMonopole-%d.ic.q'%k)
    target_mollifier(g).save('AcousticMonopole.target_mollifier.f')
    control_mollifier(g).save('AcousticMonopole.control_mollifier.f')
