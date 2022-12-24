#!/usr/bin/env python
import numpy as np
import plot3dnasa as p3d

def grid(size,alpha=1.0):
    x_min = 0.0
    x_max = alpha
    y_min = 0.0
    y_max = 1.0
    g = p3d.Grid().set_size(size, True)
    x = np.linspace(x_min, x_max, g.size[0,0] + 1)[:-1]
    y = np.linspace(y_min, y_max, g.size[0,1] + 1)[:-1]
    g.xyz[0][:,:,0,:2] = np.transpose(np.meshgrid(x, y))
    return g

def initial_condition(g, mach_number=0.1, gamma=1.4):
    y = g.xyz[0][:,:,:,1]
    x = g.xyz[0][:,:,:,0]

    temp = np.zeros_like(g.xyz[0][:,:,:,0])
    temp = 1.0 / (gamma-1)

    s = p3d.Solution().copy_from(g).quiescent(gamma)
    s.q[0][:,:,:,0] = 1.0
    s.q[0][:,:,:,1:4] = 0.0
    s.q[0][:,:,:,2] = mach_number * np.cos(2.0*np.pi*(x-0.1*np.sin(2.0*np.pi*y)))
    s.q[0][:,:,:,4] = 1.0 / gamma
    return s.fromprimitive()

def random_solution(g,time=0.0,timestep=0):

    gamma = 1.4
    s = p3d.Solution().copy_from(g).quiescent(gamma)
    n = s.get_size(0)
    s.q[0][:,:,0,0] = ( 2.0 * np.random.rand(n[0],n[1]) - 1.0 ) * 1.0e-2 + 1.0
    s.q[0][:,:,0,1] = ( 2.0 * np.random.rand(n[0],n[1]) - 1.0 ) * 1.0e-2 + 0.0
    s.q[0][:,:,0,2] = ( 2.0 * np.random.rand(n[0],n[1]) - 1.0 ) * 1.0e-2 + 0.1
    s.q[0][:,:,0,4] = ( 2.0 * np.random.rand(n[0],n[1]) - 1.0 ) * 1.0e-2 + 1./(gamma-1.)
    s.q[0][:,:,0,4] *= (gamma-1.)/gamma * s.q[0][:,:,0,0]

    s.time = time
    s._format.aux_header[0] = timestep

    return s.fromprimitive(gamma)

def target_mollifier(g, alpha = 1.0):
    x = g.xyz[0][:,:,:,0]
    y = g.xyz[0][:,:,:,1]
    x_min =  0.6  * alpha
    x_max =  1.0  * alpha
    y_min =  0.0
    y_max =  1.0
    f = p3d.Function().copy_from(g)
    f.f[0].fill(1.)
    n = f.get_size(0)
    #for i in range(n[0]):
    #    f.f[0][i,:,0,0] *= p3d.tanh_support(
    #        g.xyz[0][i,:,0,1], y_min, y_max, 40., 0.2)
    for j in range(n[1]):
        f.f[0][:,j,0,0] *= p3d.tanh_support(
            g.xyz[0][:,j,0,0], x_min, x_max, 10., 0.3)
    imin, imax = p3d.find_extents(g.xyz[0][:,0,0,0], x_min, x_max)
    jmin, jmax = p3d.find_extents(g.xyz[0][0,:,0,1], y_min, y_max)
    print ('  {:<20} {:<21} {:>4d} {:>7d}' + 6 * ' {:>4d}').format(
        'targetRegion', 'COST_TARGET', 1, 0, imin, imax, jmin, jmax, 1, -1)
    return f

def control_mollifier(g,alpha=1.0):
    x_min =  0.0  * alpha
    x_max =  0.6  * alpha
    y_min =  0.
    y_max =  1.
    f = p3d.Function().copy_from(g)
    f.f[0].fill(1.)
    n = f.get_size(0)
    for j in range(n[1]):
        f.f[0][:,j,0,0] *= p3d.tanh_support(
            g.xyz[0][:,j,0,0], x_min, x_max, 10., 0.3)
    #for i in range(n[0]):
    #    f.f[0][i,:,0,0] *= p3d.cubic_bspline_support(
    #        g.xyz[0][i,:,0,1], y_min, y_max)
    imin, imax = p3d.find_extents(g.xyz[0][:,0,0,0], x_min, x_max)
    jmin, jmax = p3d.find_extents(g.xyz[0][0,:,0,1], y_min, y_max)
    print ('  {:<20} {:<21} {:>4d} {:>7d}' + 6 * ' {:>4d}').format(
        'controlRegion', 'ACTUATOR', 1, 0, imin, imax, jmin, jmax, 1, -1)
    return f

if __name__ == '__main__':
    g = grid([256,256],alpha=0.9)
    g.save('KolmogorovFlow.xyz')
    initial_condition(g,mach_number=0.05).save('KolmogorovFlow.ic.q')
    target_mollifier(g,alpha=0.9).save('KolmogorovFlow.target_mollifier.f')
    control_mollifier(g,alpha=0.9).save('KolmogorovFlow.control_mollifier.f')
    #random_solution(g).save('ChannelFlow.ic.q')
