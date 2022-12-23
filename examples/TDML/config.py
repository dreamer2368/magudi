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

def grid(size, N=16):
    x_min =  0.
    x_max =  7.0 * N
    y_min = -3.5 * N
    y_max =  3.5 * N
    g = p3d.Grid().set_size(size, True)
    # x = x_min + 0.5 * (x_max - x_min) * (1. + mapping_function(
    #     np.linspace(0., 1., g.size[0,0]), 80., 0.58, 0.07))
    x = np.linspace(x_min, x_max, g.size[0,0] + 1)[:-1]
    y = y_min + 0.5 * (y_max - y_min) * (1. + mapping_function(
        np.linspace(0., 1., g.size[0,1]), 20., 0.62, 0.2))
    g.xyz[0][:,:,0,:2] = np.transpose(np.meshgrid(x, y))
    return g

def target_state(g, u1, u2, S, gamma=1.4):
    s = p3d.Solution().copy_from(g).quiescent(gamma)
    for i, xyz in enumerate(g.xyz):
        x = xyz[:,:,:,0]
        s.q[i][:,:,:,1] = u2 + 0.5 * (u1 - u2) * (
            1. + np.tanh(2. * xyz[:,:,:,1] / \
                         (1. + S * np.where(x > 0., x, np.zeros_like(x)))))
    return s.fromprimitive(gamma)

def initial_condition(g, u1, u2, gamma=1.4):
    s = p3d.Solution().copy_from(g).quiescent(gamma)
    for i, xyz in enumerate(g.xyz):
        s.q[i][:,:,:,1] = u2 + 0.5 * (u1 - u2) * (
            1. + np.tanh(2. * xyz[:,:,:,1]))
    return s.fromprimitive(gamma)

def target_mollifier(g):
    x_min =   0.
    x_max = 100.
    y_min = -74.
    y_max = -66.
    f = p3d.Function().copy_from(g)
    for i, xyz in enumerate(g.xyz):
        f.f[i][:,:,:,0] = p3d.tanh_support(xyz[:,:,:,0], x_min, x_max,
                                           80., 0.14) * \
            p3d.cubic_bspline_support(xyz[:,:,:,1], y_min, y_max)
        imin, imax = p3d.find_extents(xyz[:,0,0,0], x_min, x_max)
        jmin, jmax = p3d.find_extents(xyz[0,:,0,1], y_min, y_max)
        if imin and imax and jmin and jmax:
            print ('  {:<20} {:<21} {:>4d} {:>7d}' + 6 * ' {:>4d}').format(
                'targetRegion', 'COST_TARGET', i + 1, 0,
                imin, imax, jmin, jmax, 1, -1)
    return f

def control_mollifier(g):
    x_min =  1.
    x_max =  7.
    y_min = -3.
    y_max =  3.
    f = p3d.Function().copy_from(g)
    for i, xyz in enumerate(g.xyz):
        f.f[i][:,:,:,0] = p3d.tanh_support(
            xyz[:,:,:,0], x_min, x_max, 20., 0.8) * p3d.tanh_support(
                xyz[:,:,:,1], y_min, y_max, 20., 0.8)
        imin, imax = p3d.find_extents(xyz[:,0,0,0], x_min, x_max)
        jmin, jmax = p3d.find_extents(xyz[0,:,0,1], y_min, y_max)
        if imin and imax and jmin and jmax:
            print ('  {:<20} {:<21} {:>4d} {:>7d}' + 6 * ' {:>4d}').format(
                'controlRegion', 'ACTUATOR', i + 1, 0,
                imin, imax, jmin, jmax, 1, -1)
    return f

def mean_pressure(s):
    f = p3d.Function().copy_from(s)
    s = s.toprimitive()
    for i, q in enumerate(s.q):
        f.f[i][:,:,:,0] = q[:,:,:,4]
    return f

if __name__ == '__main__':
    g = grid([961, 641])
    g.save('WeiFreundSDML.xyz')
    initial_condition(g, u1=0.9, u2=0.2).save('WeiFreundSDML.ic.q')
    target_state(g, u1=0.9, u2=0.2, S=0.05).save('WeiFreundSDML.target.q')
    target_mollifier(g).save('WeiFreundSDML.target_mollifier.f')
    control_mollifier(g).save('WeiFreundSDML.control_mollifier.f')
