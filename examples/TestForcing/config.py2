#!/usr/bin/env python
import numpy as np
import plot3dnasa as p3d

def grid(size):
    x_min = 0.
    x_max = 1.
    y_min = 0.
    y_max = 1.
    g = p3d.Grid().set_size(size, True)
    x = np.linspace(x_min, x_max, g.size[0,0])
    y = np.linspace(y_min, y_max, g.size[0,1])
    g.xyz[0][:,:,0,:2] = np.transpose(np.meshgrid(x, y))
    return g

def target_mollifier(g):
    x_min =  0.
    x_max =  1.
    y_min =  0.
    y_max =  1.
    f = p3d.Function().copy_from(g)
    f.f[0].fill(1.)
    n = f.get_size(0)
    imin, imax = p3d.find_extents(g.xyz[0][:,0,0,0], x_min, x_max)
    jmin, jmax = p3d.find_extents(g.xyz[0][0,:,0,1], y_min, y_max)
    print ('  {:<20} {:<21} {:>4d} {:>7d}' + 6 * ' {:>4d}').format(
        'targetRegion', 'COST_TARGET', 1, 0, imin, imax, jmin, jmax, 1, -1)
    return f

def control_mollifier(g):
    x_min =  0.
    x_max =  1.    
    y_min =  0.
    y_max =  1.
    f = p3d.Function().copy_from(g)
    f.f[0].fill(1.)
    n = f.get_size(0)
    imin, imax = p3d.find_extents(g.xyz[0][:,0,0,0], x_min, x_max)
    jmin, jmax = p3d.find_extents(g.xyz[0][0,:,0,1], y_min, y_max)
    print ('  {:<20} {:<21} {:>4d} {:>7d}' + 6 * ' {:>4d}').format(
        'controlRegion', 'ACTUATOR', 1, 0, imin, imax, jmin, jmax, 1, -1)
    return f

def mean_pressure(s):
    f = p3d.Function().copy_from(s)
    f.f[0][:,:,:,0] = s.toprimitive().q[0][:,:,:,4]
    return f

if __name__ == '__main__':
    g = grid([4, 4])
    g.save('TestForcing.xyz')
    target_mollifier(g).save('TestForcing.target_mollifier.f')
    control_mollifier(g).save('TestForcing.control_mollifier.f')
