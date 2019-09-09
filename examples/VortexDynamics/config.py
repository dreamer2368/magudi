#!/usr/bin/env python
import numpy as np
import plot3dnasa as p3d

def grid(size, L):
    x_min = -L
    x_max =  L
    y_min = -L
    y_max =  L
    g = p3d.Grid().set_size(size, True)

    xi = np.linspace(0.,1., g.size[0,0])
    zeta = np.linspace(0.,1., g.size[0,1])

    offset1, offset2 = 0.15, 0.85
    width = 0.02

    x = 2.*L* mapping(xi, offset1, offset2, width) - L
    y = 2.*L* mapping(zeta, offset1, offset2, width) - L
    g.xyz[0][:,:,0,:2] = np.transpose(np.meshgrid(x, y))
    return g

def mapping(xi, offset1, offset2, width):
    c0 = width * np.log(np.cosh(offset2/width)/np.cosh(offset1/width))
    return xi + width * np.log( np.cosh((xi-offset2)/width)/np.cosh((xi-offset1)/width) ) - c0 + 2. * c0 * xi

def initial_condition(g, R=1.0/0.15, Ma=0.56,
                      Pr=0.7, gamma=1.4):
    from numpy.matlib import repmat
    from scipy.sparse import spdiags, kron, eye, csr_matrix
    from scipy.sparse.linalg import spsolve

    n = g.get_size(0)[0]
    # dx = ( g.xyz[0][-1:,0,0,0] - g.xyz[0][0,0,0,0] )/(n-1)
    # B2 = repmat([[1./90.], [-3./20.], [3./2.], [-49./18.], [3./2.], [-3./20.], [1./90.]],1,n)
    # D2 = spdiags(B2,[-3, -2, -1, 0, 1, 2, 3],n,n)
    # K = ( kron(eye(n),D2) + kron(D2,eye(n)) )/dx/dx
    # K = csr_matrix(K)

    xg, yg = g.xyz[0][:,:,0,0], g.xyz[0][:,:,0,1]
    src = np.zeros([n,n],dtype=np.double)
    vx = np.zeros([n,n],dtype=np.double)
    vy = np.zeros([n,n],dtype=np.double)

    loci = [[0.,R],[0.,-R],[R,0.],[-R,0.]]
    for location in loci:
        radius = np.sqrt( (xg-location[0])**2 + (yg-location[1])**2 )
        vt = 1.428/radius*( 1. - np.exp(-1.25*radius**2) )
        vt_over_r = vt/radius
        dvtdr = 3.57*np.exp(-1.25*radius**2) - 1.428/(radius**2)*( 1. - np.exp(-1.25*radius**2) )

        src += -2.0 * dvtdr * vt_over_r
        vx += -vt * (yg-location[1]) / radius
        vy += vt * (xg-location[0]) / radius

    # np.savetxt('loci.txt', loci)
    #
    # import subprocess
    # command = 'matlab -nodesktop -nosplash -r "poisson(' + str(n) + ',' + str(R) + ',' + str(Ma) + ',' + str(gamma) + '); exit"'
    # subprocess.check_call(command,shell=True)
    # p0 = np.fromfile('p0.bin',dtype=np.double)
    # rho0 = np.fromfile('rho0.bin',dtype=np.double)
    # p0, rho0 = np.reshape(p0,[n,n]).T, np.reshape(rho0,[n,n]).T

    s = p3d.Solution().copy_from(g).quiescent(gamma)
    for i, xyz in enumerate(g.xyz):
        # s.q[i][:,:,0,0] = rho0
        s.q[i][:,:,0,0] = 1.0
        s.q[i][:,:,0,1] = vx
        s.q[i][:,:,0,2] = vy
        # s.q[i][:,:,0,4] = p0
        s.q[i][:,:,0,4] = 1.0/gamma
    return s.fromprimitive(gamma)

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

if __name__ == '__main__':
    R = 1.0/0.15
    g = grid([301, 301], 10.*R)
    g.save('VortexDynamics.xyz')
    initial_condition(g).save('VortexDynamics.ic.q')
    # target_mollifier(g).save('AcousticMonopole.target_mollifier.f')
    # control_mollifier(g).save('AcousticMonopole.control_mollifier.f')
