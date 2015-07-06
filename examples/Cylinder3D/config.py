#!/usr/bin/env python
from pylab import *
import numpy as np
import plot3dnasa as p3d

def mapping_function(x, sigma):
    return np.sinh(sigma * x) / np.sinh(sigma)

def get_mesh(x1, x2, a, b, xc, sigma):
    f = lambda x: a + b * np.tanh((x - xc) / sigma)
    x = np.array([x1])
    while x[-1] < x2:
        x = np.append(x, [x[-1] + f(x[-1])])
    if np.abs(x[-1] + f(x[-1]) - x2) < np.abs(x[-1] - x2):
        x = np.append(x, [x[-1] + f(x[-1])])        
    x = x1 + (x2 - x1) * (x - x1) / (x[-1] - x1)
    return x

def max_stretching(sigma, x1, x2, dx1, dx2):
    assert x2 > x1
    xc = (x2 * dx1 + x1 * dx2) / (dx1 + dx2)
    b = (dx2 - dx1) / (np.tanh((x2 - xc) / sigma) -
                       np.tanh((x1 - xc) / sigma))
    a = 0.5 * (dx1 + dx2 - b * (np.tanh((x1 - xc) / sigma) +
                                np.tanh((x2 - xc) / sigma)))
    x = get_mesh(x1, x2, a, b, xc, sigma)
    return np.amax(np.abs(x[2:] - 2. * x[1:-1] + x[:-2]) / (x[1:-1] - x[:-2]))

def spacing_func(x, dx1, dx2, sigma):
    x1 = x.min()
    x2 = x.max()
    xc = (x2 * dx1 + x1 * dx2) / (dx1 + dx2)
    b = (dx2 - dx1) / (np.tanh((x2 - xc) / sigma) -
                       np.tanh((x1 - xc) / sigma))
    a = 0.5 * (dx1 + dx2 - b * (np.tanh((x1 - xc) / sigma) +
                                np.tanh((x2 - xc) / sigma)))
    return a + b * np.tanh((x - xc) / sigma)

def mesh_segment(x1, x2, dx1, dx2, s_max):
    from scipy.optimize import fsolve
    sigma = fsolve(lambda x: (max_stretching(x, x1, x2, dx1, dx2) -
                              s_max) ** 2, 2.)
    xc = (x2 * dx1 + x1 * dx2) / (dx1 + dx2)
    b = (dx2 - dx1) / (np.tanh((x2 - xc) / sigma) -
                       np.tanh((x1 - xc) / sigma))
    a = 0.5 * (dx1 + dx2 - b * (np.tanh((x1 - xc) / sigma) +
                                np.tanh((x2 - xc) / sigma)))
    return get_mesh(x1, x2, a, b, xc, sigma)

def radial_coordinate_Inoue(r_breakpts=[0.5, 0.95, 20., 95., 350., 800.],
                      dr=[0.005, 0.2, 20.]):
    r = np.linspace(r_breakpts[0], r_breakpts[1], int(np.rint(
        (r_breakpts[1] - r_breakpts[0]) / dr[0])) + 1)[:-1]
    r = np.append(r, mesh_segment(r_breakpts[1], r_breakpts[2],
                                  dr[0], dr[1], 0.04)[:-1])
    r = np.append(r, np.linspace(r_breakpts[2], r_breakpts[3], int(np.rint(
        r_breakpts[3] - r_breakpts[2]) / dr[1]) + 1)[:-1])
    r = np.append(r, mesh_segment(r_breakpts[3], r_breakpts[4],
                                  dr[1], dr[2], 0.09)[:-1])
    r = np.append(r, np.linspace(r_breakpts[4], r_breakpts[5], int(np.rint(
        r_breakpts[5] - r_breakpts[4]) / dr[2]) + 1))
    return r

def radial_coordinate(n, mapping_type='geom', r_buffer=800., dr_min=0.005,
                      num_uniform=50):
    from scipy.optimize import fsolve
    r_min = 0.5
    r_max = r_buffer
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

def draw_line3d(f, label, p1, p2):
    print >>f, 'gg::conBegin'
    print >>f, 'gg::segBegin -type 3D_LINE'
    print >>f, 'gg::segAddControlPt {%f %f %f}' % (p1[0], p1[1], p1[2])
    print >>f, 'gg::segAddControlPt {%f %f %f}' % (p2[0], p2[1], p2[2])
    print >>f, 'gg::segEnd'
    print >>f, 'set %s [gg::conEnd]' % (label)
    return None

def write_glyph_file(filename, grid_size, r_surface=1., r_sound=100.,
                     r_buffer=1500., n_surface=70, n_buffer=340,
                     dr_min=0.005, dr_sound=0.2,
                     dr_buffer=20.):
    # Write glyph file.
    with open(filename, 'w') as f:
        print >>f, """
        gg::memClear
        gg::defReset
        gg::tolReset
        gg::dispViewReset
        """
        draw_line3d(f, 'tmp_wake', [0.5, 0., 0.], [r_buffer, 0., 0.])
        print >>f, 'gg::conDim [list $tmp_wake] %i' % (grid_size[0])
        print >>f, 'gg::conSetBreakPt $tmp_wake {%f 0 0}' % r_surface
        print >>f, 'gg::conSetBreakPt $tmp_wake {%f 0 0}' % r_sound
        print >>f, 'gg::conSubConDim $tmp_wake [list %i %i %i]' % \
            (n_surface + 1, grid_size[0] - (n_surface + n_buffer),
             n_buffer + 1)
        print >>f, 'gg::conDistFunc $tmp_wake -function MRQS'
        print >>f, 'gg::conBeginSpacing $tmp_wake -sub 1 %f' % dr_min
        print >>f, 'set ds_end [ggu::vec3Length [ggu::vec3Sub ' \
            '[gg::conGetPt $tmp_wake %i] ' \
            '[gg::conGetPt $tmp_wake %i]]]' % (n_surface, n_surface + 1)
        print >>f, 'gg::conEndSpacing $tmp_wake -sub 1 $ds_end'
        print >>f, 'gg::conBeginSpacing $tmp_wake -sub 2 $ds_end'
        print >>f, 'gg::conEndSpacing $tmp_wake -sub 2 %f' % dr_sound
        print >>f, 'gg::conBeginSpacing $tmp_wake -sub 3 %f' % dr_sound
        print >>f, 'gg::conEndSpacing $tmp_wake -sub 3 %f' % dr_buffer
        print >>f, """
        gg::domExtrusionBegin $tmp_wake -default ROTATE
        gg::domExtrusionAtt -angle 360
        gg::domExtrusionAtt -axis {0 0 0} {0 0 1}
        """
        print >>f, 'gg::domExtrusionStep %i' % (grid_size[1] - 1)
        print >>f, 'set dom_cylinder [gg::domExtrusionEnd]'
        print >>f, 'gg::domTransformBegin $dom_cylinder -maintain_linkage'
        print >>f, 'gg::xformRotate [list 0 0 0] [list 0 0 1] %f' % \
            (90. + 180. / grid_size[1])
        print >>f, 'gg::domTransformEnd'
        print >>f, 'gg::domExport $dom_cylinder \"Cylinder.xyz\" ' \
            '-style PLOT3D -form UNFORMATTED -precision DOUBLE -endian NATIVE'

def grid_from_gridgen(**kwargs):
    if not os.path.isfile('Cylinder.xyz'):
        f = tempfile.NamedTemporaryFile(delete=False)
        write_glyph_file(f.name, **kwargs)
        subprocess.check_output(["gridgen", "-b", f.name])
        os.unlink(f.name)

def grid(num_polar, num_axial=1):
    # r = radial_coordinate_Inoue()
    r = radial_coordinate(299, dr_min=0.0005, r_buffer=80., num_uniform=12)
    theta = np.linspace(-np.pi, np.pi, num_polar)
    z = np.linspace(0., 2., num_axial + 1)[:-1]
    g = p3d.Grid().set_size([r.size, theta.size, z.size], True)
    r, theta = np.meshgrid(r, theta)
    for k in range(num_axial):
        g.xyz[0][:,:,k,0] = np.transpose(r * np.cos(theta))
        g.xyz[0][:,:,k,1] = np.transpose(r * np.sin(theta))
        g.xyz[0][:,:,k,2] = z[k]
    return g

def target_state(g, mach_number, gamma=1.4):
    s = p3d.Solution().copy_from(g).quiescent(gamma)
    s.q[0][:,:,:,1] = mach_number
    return s.fromprimitive(gamma)

def target_mollifier(g):
    r_min =  40.
    r_max =  60.
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
    x_min = 0.7
    x_max = 1.3
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

if __name__ == '__main__':
    g = grid(num_polar=501, num_axial=96)
    g.save('Cylinder.xyz')
    target_state(g, mach_number=0.2).save('Cylinder.target.q')
    target_mollifier(g).save('Cylinder.target_mollifier.f')
    control_mollifier(g).save('Cylinder.control_mollifier.f')
