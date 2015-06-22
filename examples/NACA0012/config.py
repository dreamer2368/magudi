#!/usr/bin/env python
import numpy as np
import plot3dnasa as p3d
import tempfile
import subprocess
import os
import os.path

def airfoil_profile(x, c, t):
    return 5. * t * c * (0.2969 * np.sqrt(x / c) - 0.1260 * (x / c) -
                         0.3516 * (x / c) ** 2 + 0.2843 * (x / c) ** 3 -
                         0.1015 * (x / c) ** 4)

def airfoil_slope(x, c, t):
    return 5. * t * (0.14845 / np.sqrt(x / c) - 0.1260 - 0.7032 * (x / c) +
                     0.8529 * (x / c) ** 2 - 0.406 * (x / c) ** 3)

def write_glyph_file(filename, grid_size, t=0.12, alpha=2., ds_te=0.0001,
                     ds_le=0.0008, ds_normal=0.00025, stretch_ratio=1.02,
                     stop_height=100.):
    c = 2000.
    s = np.linspace(0., 1., 200)
    x_airfoil = np.cumsum(np.tanh(3.5 * (1. - s)) + np.tanh(3.5 * s) - 1.)
    x_airfoil = c * (x_airfoil - x_airfoil.min()) / \
                (x_airfoil.max() - x_airfoil.min())
    y_airfoil = airfoil_profile(x_airfoil, c, t)
    x_te = x_airfoil[-1] - 0.5 * y_airfoil[-1] / \
           airfoil_slope(x_airfoil[-1], c, t)
    x_start = 0.52 * c
    # Write glyph file.
    with open(filename, 'w') as f:
        print >>f, """
        gg::memClear
        gg::defReset
        gg::tolReset
        gg::dispViewReset
        """
        print >>f, 'gg::dbCurveBegin -type CUBIC'
        for x, y in zip(x_airfoil[1:], y_airfoil[1:]):
            print >>f, 'gg::dbCurveAddPt {%.6f %.6f %.6f}' % (x, y, 0.)
        print >>f, 'set db_NACA1 [gg::dbCurveEnd]'
        print >>f, 'gg::dbCurveBegin -type CUBIC'
        for x, y in zip(reversed(x_airfoil[1:]), reversed(y_airfoil[1:])):
            print >>f, 'gg::dbCurveAddPt {%.6f %.6f %.6f}' % (x, -y, 0.)
        print >>f, 'set db_NACA2 [gg::dbCurveEnd]'
        print >>f, 'gg::dbCurveBegin -type CONIC -rho 0.5'
        print >>f, 'gg::dbCurveAddPt {%.6f %.6f %.6f}' % \
            (x_airfoil[1],  y_airfoil[1], 0.)
        print >>f, 'gg::dbCurveAddPt {%.6f %.6f %.6f}' % \
            (x_airfoil[1], -y_airfoil[1], 0.)
        print >>f, 'gg::dbCurveAddPt {%.6f %.6f %.6f}' \
            % (x_airfoil[0],  y_airfoil[0], 0.)
        print >>f, 'set db_NACA3 [gg::dbCurveEnd]'
        print >>f, 'gg::dbCurveBegin -type CONIC -rho 0.5'
        print >>f, 'gg::dbCurveAddPt {%.6f %.6f %.6f}' % \
            (x_airfoil[-1],  y_airfoil[-1], 0.)
        print >>f, 'gg::dbCurveAddPt {%.6f %.6f %.6f}' % \
            (x_airfoil[-1], -y_airfoil[-1], 0.)
        print >>f, 'gg::dbCurveAddPt {%.6f %.6f %.6f}' % \
            (x_te,  0., 0.)
        print >>f, 'set db_NACA4 [gg::dbCurveEnd]'
        print >>f, """
        set con_NACA1 [gg::conOnDBEnt $db_NACA1]
        set con_NACA2 [gg::conOnDBEnt $db_NACA2]
        set con_NACA3 [gg::conOnDBEnt $db_NACA3]
        set con_NACA4 [gg::conOnDBEnt $db_NACA4]
        """
        print >>f, 'set con_NACA5 [gg::conSplit $con_NACA1 ' \
            '[gg::conGetPt $con_NACA1 -x %f]]' % (x_start)
        print >>f, """
        set con_NACA6 [gg::conJoin $con_NACA5 $con_NACA4]
        set con_NACA7 [gg::conJoin $con_NACA6 $con_NACA2]
        set con_NACA8 [gg::conJoin $con_NACA7 $con_NACA3]
        set con_NACA [gg::conJoin $con_NACA1 $con_NACA8]
        """
        print >>f, 'gg::conDim $con_NACA %i' % (grid_size[0])
        print >>f, 'set LE_pt [gg::conGetPt $con_NACA -x 0]'
        print >>f, 'gg::conSetBreakPt $con_NACA $LE_pt'
        print >>f, 'set TE_pt [gg::conGetPt $con_NACA -y 0]'
        print >>f, 'gg::conSetBreakPt $con_NACA $TE_pt'
        print >>f, 'gg::conBreakPtSpacing $con_NACA -breakpt 1 %f' % \
            (ds_te * c)
        print >>f, 'gg::conBreakPtSpacing $con_NACA -breakpt 2 %f' % \
            (ds_le * c)
        print >>f, 'set ds_start [ggu::vec3Length [ggu::vec3Sub [' \
            'gg::conGetPt $con_NACA %i] [gg::conGetPt $con_NACA %i]]]' % \
            (grid_size[0], grid_size[0] - 1)
        print >>f, 'gg::conBeginSpacing $con_NACA -sub 1 $ds_start'
        print >>f, 'gg::conEndSpacing $con_NACA -sub 3 $ds_start'
        print >>f, 'gg::domExtrusionBegin [list $con_NACA] -edge'
        print >>f, 'gg::domExtrusionMode HYPERBOLIC'
        print >>f, 'gg::domExtrusionAtt -s_init %f -march_plane {0 0 1} ' \
            '-stop_height %f -growth_geometric %f -normal_count 20 ' \
            '-normal_relax 1 -vol_smoothing 0.25' % \
            (ds_normal * c, stop_height * c, stretch_ratio)
        print >>f, 'gg::domExtrusionStep -result ExtResult %i' % \
            (grid_size[1] - 1)
        print >>f, 'set dom_NACA [gg::domExtrusionEnd]'
        print >>f, 'gg::domTransformBegin $dom_NACA -maintain_linkage'
        print >>f, 'gg::xformRotate [list 0 0 0] [list 0 0 1] %f' % (-alpha)
        print >>f, 'gg::domTransformEnd'
        print >>f, 'gg::domTransformBegin $dom_NACA -maintain_linkage'
        print >>f, 'gg::xformScale [list 0 0 0] [list %f %f %f]' % \
            (1. / x_te, 1. / x_te, 1. / x_te)
        print >>f, 'gg::domTransformEnd'
        print >>f, 'gg::dbDelete ALL'
        print >>f, 'gg::domExport $dom_NACA \"NACA0012.xyz\" -style PLOT3D ' \
            '-form UNFORMATTED -precision DOUBLE -endian NATIVE'

def grid(**kwargs):
    if not os.path.isfile('NACA0012.xyz'):
        f = tempfile.NamedTemporaryFile(delete=False)
        write_glyph_file(f.name, **kwargs)
        subprocess.check_output(["gridgen", "-b", f.name])
        os.unlink(f.name)

def target_state(g, mach_number, gamma=1.4):
    s = p3d.Solution().copy_from(g).quiescent(gamma)
    s.q[0][:,:,:,1] = mach_number
    return s.fromprimitive(gamma)

if __name__ == '__main__':
    grid(grid_size = [512, 448])
    g = p3d.fromfile('NACA0012.xyz')
    target_state(g, mach_number=0.5).save('NACA0012.target.q')
