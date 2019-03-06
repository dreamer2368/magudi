#!/usr/bin/env python
import numpy as np
from scipy.optimize import fsolve
import numpy.random as random

import plot3dnasa as p3d

EPSILON = np.finfo(float).eps

def grid(numX,numY):
    x_min = -14.
    x_max =  14.
    y_min = -14.
    y_max =  14.
    size = [[numX,numY],[numX,numY]]
    g = p3d.Grid().set_size(size,True)
    x1 = np.linspace(x_min, 0., g.size[0,0])
    x2 = np.linspace(0., x_max, g.size[1,0])
    y = np.linspace(y_min, y_max, g.size[0,1])
    g.xyz[0][:,:,0,2] = 0.
    g.xyz[1][:,:,0,2] = 0.
    g.xyz[0][:,:,0,:2] = np.transpose(np.meshgrid(x1, y))
    g.xyz[1][:,:,0,:2] = np.transpose(np.meshgrid(x2, y))
    return g

def initial_condition(g, mach_number=1.3, gamma=1.4):
    s = p3d.Solution().copy_from(g).quiescent(gamma)
    return s.fromprimitive(gamma)

def target_mollifier(g):
    x_min =  -3.
    x_max =   3.
    y_min =  -1.
    y_max =   1.
    f = p3d.Function().copy_from(g)
    y = g.xyz[0][0,:,0,1]
    n = f.get_size()
    block_code = ['E','W']
    for i, fi in enumerate(f.f):
        x = g.xyz[i][:,0,0,0]
        if x.max() < x_min or x.min() > x_max:
            fi.fill(0.)
            continue
        fi.fill(1.)
        for j in range(n[i][0]):
            fi[j,:,0,0] *= p3d.tanh_support(y, y_min, y_max, 5., 0.2)
        for k in range(n[i][1]):
            fi[:,k,0,0] *= p3d.cubic_bspline_support(x, x_min, x_max,strict=False)
        imin, imax = p3d.find_extents(x, x_min, x_max)
        jmin, jmax = p3d.find_extents(y, y_min, y_max)
        print ('  {:<20} {:<21} {:>4d} {:>7d}' + 6 * ' {:>4d}').format(
            'targetRegion'+block_code[i], 'COST_TARGET', i+1, 0, imin, imax, jmin, jmax, 1, -1)
    return f

def control_mollifier(g):
    x_min = -3.
    x_max =  3.    
    y_min =  2.
    y_max =  4.
    f = p3d.Function().copy_from(g)
    y = g.xyz[0][0,:,0,1]
    n = f.get_size()
    block_code = ['E','W']
    for i, fi in enumerate(f.f):
        x = g.xyz[i][:,0,0,0]
        if x.max() < x_min or x.min() > x_max:
            fi.fill(0.)
            continue
        fi.fill(1.)
        for j in range(n[i][0]):
            fi[j,:,0,0] *= p3d.cubic_bspline_support(y, y_min, y_max)
        for k in range(n[i][1]):
            fi[:,k,0,0] *= p3d.cubic_bspline_support(x, x_min, x_max,strict=False)
        imin, imax = p3d.find_extents(x, x_min, x_max)
        jmin, jmax = p3d.find_extents(y, y_min, y_max)
        print ('  {:<20} {:<21} {:>4d} {:>7d}' + 6 * ' {:>4d}').format(
            'controlRegion'+block_code[i], 'ACTUATOR', i+1, 0, imin, imax, jmin, jmax, 1, -1)
    return f

def mean_pressure(s):
    f = p3d.Function().copy_from(s)
    for i, fi in enumerate(f.f):
        fi[:,:,:,0] = s.toprimitive().q[i][:,:,:,4]
    return f


if __name__ == '__main__':
    g = grid(101,201)
    g.save('MultiblockAcousticMonopole.xyz')
    initial_condition(g).save('MultiblockAcousticMonopole.ic.q')
    control_mollifier(g).save('MultiblockAcousticMonopole.control_mollifier.f')
    target_mollifier(g).save('MultiblockAcousticMonopole.target_mollifier.f')
#    target_state(g).save('MultiblockJet.target.q')
#    initial_condition(g).save('MultiblockJet.ic.q')
#    gi = extract_inflow(g)
#    gi.save('MultiblockJet.inflow.xyz')
#    modes = eigenmodes()
#    for i, mode in enumerate(modes):
#        sr, si = inflow_perturbations(gi, mode)
#        sr.save('MultiblockJet-%02d.eigenmode_real.q' % (i + 1))
#        si.save('MultiblockJet-%02d.eigenmode_imag.q' % (i + 1))
