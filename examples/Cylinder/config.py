#!/usr/bin/env python
import numpy as np
import plot3dnasa as p3d

def mapping_function(x, sigma):
    return np.sinh(sigma * x) / np.sinh(sigma)

def radial_coordinate(n, mapping_type='geom', r_buffer=800., dr_min=0.005):
    from scipy.optimize import fsolve
    r_min = 0.5
    r_max = r_buffer
    num_uniform = 50
    if mapping_type == 'sinh':
        sigma = fsolve(lambda x: (
                (r_max - r_min - num_uniform * dr_min) * mapping_function(
                    1. / (n - 1.), x) - dr_min) ** 2, 2.)
        r = np.append(np.linspace(r_min, r_min + dr_min * num_uniform,
                                  num_uniform + 1), r_min + dr_min *
                      num_uniform + (r_max - r_min - num_uniform * dr_min) *
                      mapping_function(np.linspace(
                          0., 1., n - num_uniform), sigma)[1:])
    else:
        sigma = fsolve(lambda x: (r_max - r_min) / dr_min - num_uniform + 1 -
                       (x ** (n - num_uniform) - 1.) /
                       (x - 1.), 1.02)
        r = r_min + np.append([0.], np.cumsum(
            [dr_min if r < num_uniform - 1
             else dr_min * sigma ** (r - num_uniform + 1)
             for r in range(n - 1)]))
    return r    

def grid(grid_size):
    r = radial_coordinate(grid_size[0], r_buffer=300.)
    print r[-57]
    theta = np.linspace(-np.pi, np.pi, grid_size[1])
    g = p3d.Grid().set_size([r.size, theta.size, 1], True)
    r, theta = np.meshgrid(r, theta)
    g.xyz[0][:,:,0,0] = np.transpose(r * np.cos(theta))
    g.xyz[0][:,:,0,1] = np.transpose(r * np.sin(theta))
    return g

def target_state(g, mach_number, gamma=1.4):
    s = p3d.Solution().copy_from(g).quiescent(gamma)
    s.q[0][:,:,:,1] = mach_number
    return s.fromprimitive(gamma)

def target_mollifier(g):
    r_min =  60.
    r_max =  90.
    r = np.sqrt(g.xyz[0][:,0,0,0] ** 2 + g.xyz[0][:,0,0,1] ** 2)
    f = p3d.Function().copy_from(g)
    f.f[0].fill(1.)
    n = f.get_size(0)
    for j in range(n[1]):
        f.f[0][:,j,0,0] *= p3d.cubic_bspline_support(r, r_min, r_max)
    imin, imax = p3d.find_extents(r, r_min, r_max)
    print ('  {:<20} {:<21} {:>4d} {:>7d}' + 6 * ' {:>4d}').format(
        'targetRegion', 'COST_TARGET', 1, 0, imin, imax, 1, -1, 1, -1)
    return f

def control_mollifier(g):
    x_min = 0.8
    x_max = 1.2
    y_min = -0.25
    y_max = 0.25
    r_min =  x_min
    r_max =  np.sqrt(x_max ** 2 + y_max ** 2)
    theta_min = np.arctan2(y_min, x_min)
    theta_max = np.arctan2(y_max, x_min)
    r = np.sqrt(g.xyz[0][:,:,0,0] ** 2 + g.xyz[0][:,:,0,1] ** 2)
    theta = np.arctan2(g.xyz[0][:,:,0,1], g.xyz[0][:,:,0,0])
    f = p3d.Function().copy_from(g)
    f.f[0].fill(1.)
    n = f.get_size(0)
    f.f[0][:,:,0,0] *= p3d.tanh_support(g.xyz[0][:,:,0,0], x_min, x_max,
                                        10., 0.3)
    f.f[0][:,:,0,0] *= p3d.tanh_support(g.xyz[0][:,:,0,1], y_min, y_max,
                                        10., 0.3)
    imin, imax = p3d.find_extents(r[:,0], r_min, r_max)
    jmin = np.argmin(np.abs(theta[0,:] - theta_min))
    jmax = np.argmin(np.abs(theta[0,:] - theta_max)) + 2
    print ('  {:<20} {:<21} {:>4d} {:>7d}' + 6 * ' {:>4d}').format(
        'controlRegion', 'ACTUATOR', 1, 0, imin, imax, jmin, jmax, 1, -1)
    return f

def mean_pressure(s):
    f = p3d.Function().copy_from(s)
    f.f[0][:,:,:,0] = s.toprimitive().q[0][:,:,:,4]
    return f

if __name__ == '__main__':
    g = grid([408, 501])
    g.save('Cylinder.xyz')
    target_state(g, mach_number=0.2).save('Cylinder.target.q')
    target_mollifier(g).save('Cylinder.target_mollifier.f')
    control_mollifier(g).save('Cylinder.control_mollifier.f')
