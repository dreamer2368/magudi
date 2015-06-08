#ini!/usr/bin/env python
import numpy as np
import plot3dnasa as p3d

def mapping_function(s, b, c, sigma):
    from scipy.special import erf
    return ((s - 0.5) * (1.0 + 2.0 * b) - b * sigma *
            (np.exp(- ((s - 0.5 + c) / sigma) ** 2) / np.sqrt(np.pi) +
             ((s - 0.5 + c) / sigma) * erf((s - 0.5 + c) / sigma) -
             np.exp(- ((s - 0.5 - c) / sigma) ** 2) / np.sqrt(np.pi) -
             ((s - 0.5 - c) / sigma) * erf((s - 0.5 - c) / sigma))) / \
        (0.5 + b - b * sigma * (np.exp(- ((0.5 + c) / sigma) ** 2) /
                                np.sqrt(np.pi) + ((0.5 + c) / sigma) *
                                erf((0.5 + c) / sigma) -
                                np.exp(- ((0.5 - c) / sigma) ** 2) /
                                np.sqrt(np.pi) - ((0.5 - c) / sigma) *
                                erf((0.5 - c) / sigma)))

def initial_condition(g, u1, u2, gamma=1.4):
    s = p3d.Solution().copy_from(g).quiescent(gamma)
    s.q[0][:,:,:,1] = u2 + 0.5 * (u1 - u2) * (
        1. + np.tanh(2. * g.xyz[0][:,:,:,1]))
    return s.fromprimitive(gamma)

def target_state(g, u1, u2, S, gamma=1.4):
    s = p3d.Solution().copy_from(g).quiescent(gamma)
    x = g.xyz[0][:,:,:,0]
    s.q[0][:,:,:,1] = u2 + 0.5 * (u1 - u2) * (
        1. + np.tanh(2. * g.xyz[0][:,:,:,1] / \
                     (1. + S * np.where(x > 0., x, np.zeros_like(x)))))
    return s.fromprimitive(gamma)

def grid(size):
    x_min =  -60.
    x_max =  160.
    y_min = -100.
    y_max =  100.
    g = p3d.Grid().set_size(size, True)
    x = x_min + 0.5 * (x_max - x_min) * (1. + mapping_function(
        np.linspace(0., 1., g.size[0,0]), 80., 0.58, 0.07))
    y = y_min + 0.5 * (y_max - y_min) * (1. + mapping_function(
        np.linspace(0., 1., g.size[0,1]), 20., 0.62, 0.2))
    g.xyz[0][:,:,0,:2] = np.transpose(np.meshgrid(x, y))
    return g

def target_mollifier(g):
    x_min =   0.
    x_max = 100.
    y_min = -80.
    y_max = -60.
    f = p3d.Function().copy_from(g)
    f.f[0].fill(1.)
    n = f.get_size(0)
    for j in range(n[1]):
        f.f[0][:,j,0,0] *= p3d.tanh_support(
            g.xyz[0][:,j,0,0], x_min, x_max, 40., 0.2)    
    for i in range(n[0]):
        f.f[0][i,:,0,0] *= p3d.cubic_bspline_support(
            g.xyz[0][i,:,0,1], y_min, y_max)
    imin, imax = p3d.find_extents(g.xyz[0][:,0,0,0], x_min, x_max)        
    jmin, jmax = p3d.find_extents(g.xyz[0][0,:,0,1], y_min, y_max)
    print ('  {:<20} {:<21} {:>4d} {:>7d}' + 6 * ' {:>4d}').format(
        'targetRegion', 'COST_TARGET', 1, 0, imin, imax, jmin, jmax, 1, -1)
    return f

def control_mollifier(g):
    x_min =  1.
    x_max =  7.    
    y_min = -3.
    y_max =  3.
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

if __name__ == '__main__':
    g = grid([960, 640])
    g.save('WeiFreundSDML.xyz')
    initial_condition(g, u1=0.9, u2=0.2).save('WeiFreundSDML.ic.q')
    target_state(g, u1=0.9, u2=0.2, S=0.05).save('WeiFreundSDML.target.q')
    target_mollifier(g).save('WeiFreundSDML.target_mollifier.f')
    control_mollifier(g).save('WeiFreundSDML.control_mollifier.f')
