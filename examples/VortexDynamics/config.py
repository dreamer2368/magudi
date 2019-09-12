#!/usr/bin/env python
import numpy as np
import plot3dnasa as p3d

def grid(Nx, L, grid_ratio, length_ratio, offset_ratio):
    # Nx: number of grid points in one direction
    # L: half length of domain in one direction
    # grid_ratio: ratio of grid sizes in far-field(dx1) to near-field(dx2)
    # length_ratio: ratio of domain length of near-field(L2) to 2*L
    # offset_ratio: ratio of offset length to 2*L

    a, b, c = grid_ratio, offset_ratio, length_ratio

    dx2 = ( 1. + 2.*b*(a - 1.)/(a + 1.) + c*(a - 1.) ) / a
    dx1 = a * dx2
    dxi = 2. * b / (dx1 + dx2)
    xin = c / dx2

    inc = 1./(Nx-1)
    xi = ( np.arange(Nx-1) + 0.5 ) * inc
    xi1 = 0.5 - 0.5 * xin - dxi
    xi2 = 0.5 - 0.5 * xin
    xi3 = 0.5 + 0.5 * xin
    xi4 = 0.5 + 0.5 * xin + dxi
    cond1 = np.where(xi<=xi1)
    cond2 = np.where((xi>=xi1) & (xi<=xi2) == True)
    cond3 = np.where((xi>=xi2) & (xi<=xi3) == True)
    cond4 = np.where((xi>=xi3) & (xi<=xi4) == True)
    cond5 = np.where(xi>=xi4)

    dx = np.zeros_like(xi)
    dx[cond1] = dx1
    dx[cond2] = dx1 - (dx1-dx2) / dxi * ( xi[cond2] - 0.5 + 0.5 * xin + dxi )
    dx[cond3] = dx2
    dx[cond4] = dx1 + (dx1-dx2) / dxi * ( xi[cond4] - 0.5 - 0.5 * xin - dxi )
    dx[cond5] = dx1
    dx *= inc

    x = np.zeros((Nx,),dtype=np.double)
    for k, dxk in enumerate(dx):
        x[k+1] = x[k] + dxk

    x /= x[-1]
    x = L * ( 2. * x - 1. )
    y = x
    g = p3d.Grid().set_size([Nx,Nx], True)

    # N2 = int(np.ceil( grid_ratio/( grid_ratio - 1. + 1./length_ratio )*Nx ))
    # N1 = Nx - N2
    # dx1, dx2 = (1.-length_ratio)*L/N1, length_ratio*L/N2
    # x = np.zeros([2*(N1+N2)+1,],dtype=np.double)
    # x[:N1] = - L + dx1 * np.arange(N1)
    # x[-N1:] = L - dx1 * np.arange(N1)
    # x[-(N1+N2):-N1] = dx2 + dx2 * np.arange(N2)
    # x[N1:(N1+N2)] = - dx2 * N2 + dx2 * np.arange(N2)
    # y = x
    # g = p3d.Grid().set_size([2*(N1+N2)+1,2*(N1+N2)+1], True)

    # xi = np.linspace(0.,1., g.size[0,0])
    # zeta = np.linspace(0.,1., g.size[0,1])
    #
    # offset1, offset2 = 0.15, 0.85
    # width = 0.02
    #
    # x = 2.*L* mapping(xi, offset1, offset2, width) - L
    # y = 2.*L* mapping(zeta, offset1, offset2, width) - L
    g.xyz[0][:,:,0,:2] = np.transpose(np.meshgrid(x, y))
    return g

def mapping(xi, offset1, offset2, width):
    c0 = width * np.log(np.cosh(offset2/width)/np.cosh(offset1/width))
    return xi + 2.39 * ( width * np.log( np.cosh((xi-offset2)/width)/np.cosh((xi-offset1)/width) ) - c0 + 2. * c0 * xi )

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

    eps = 0.0e-16
    loci = [[0.,0.]]
    for location in loci:
        radius = np.sqrt( (xg-location[0])**2 + (yg-location[1])**2 )
        vt = 1.428 * Ma / (radius+eps) * ( 1. - np.exp(-1.25*radius**2) )

#        vt_over_r = vt/radius
#        dvtdr = 3.57*np.exp(-1.25*radius**2) - 1.428/(radius**2)*( 1. - np.exp(-1.25*radius**2) )
#        src += -2.0 * dvtdr * vt_over_r

        vx += -vt * (yg-location[1]) / (radius+eps)
        vy += vt * (xg-location[0]) / (radius+eps)

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
    g = grid(429, 115.*R, 100., 0.03, 0.2)
    # g = grid([429, 429], 115.*R)
    g.save('VortexDynamics.xyz')
    initial_condition(g).save('VortexDynamics.ic.q')
    # target_mollifier(g).save('AcousticMonopole.target_mollifier.f')
    # control_mollifier(g).save('AcousticMonopole.control_mollifier.f')
