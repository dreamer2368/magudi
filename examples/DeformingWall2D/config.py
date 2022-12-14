#!/usr/bin/env python
import numpy as np
import scipy.sparse as sps
import plot3dnasa as p3d

def mapping_function(x, sigma):
    return np.sinh(sigma * x) / np.sinh(sigma)

def grid(size, mapping_type='sinh'):
    from scipy.optimize import fsolve
    x_min =  -10.
    x_max =  100.
    y_min =    0.
    y_max =   25.
    z_min =  -40.
    z_max =   40.
    dy_min = 0.016
    num_uniform = 7
    g = p3d.Grid().set_size(size, True)
    x = np.linspace(x_min, x_max, g.size[0,0] + 1)[:-1]
    if mapping_type == 'sinh':
        sigma = fsolve(lambda x: (
                (y_max - y_min - num_uniform * dy_min) * mapping_function(
                    1. / (g.size[0,1] - 1.), x) - dy_min) ** 2, 2.)
        y = np.append(np.linspace(y_min, y_min + dy_min * num_uniform,
                                  num_uniform + 1), y_min + dy_min *
                      num_uniform + (y_max - y_min - num_uniform * dy_min) *
                      mapping_function(np.linspace(0., 1., g.size[0,1] -
                                                   num_uniform), sigma)[1:])
    else:
        sigma = fsolve(lambda x: (y_max - y_min) / dy_min - num_uniform + 1 -
                       (x ** (g.size[0,1] - num_uniform) - 1.) /
                       (x - 1.), 1.02)
        print 100. * (sigma - 1.)
        y = np.append([0.], np.cumsum(
                [dy_min if r < num_uniform - 1
                 else dy_min * sigma ** (r - num_uniform + 1)
                 for r in range(g.size[0,1] - 1)]))
    # z = np.linspace(z_min, z_max, g.size[0,2] + 1)[:-1]
    for i in range(x.size):
        g.xyz[0][i,:,:,0] = x[i]
    for j in range(y.size):
        g.xyz[0][:,j,:,1] = y[j]
    # for k in range(z.size):
    #     g.xyz[0][:,:,k,2] = z[k]
    return g

def solve_blasius(eta_max, Ng = 1001):
    eta = np.linspace(0., eta_max, Ng)
    dx = eta_max / Ng
    f = np.zeros(Ng,)

    # Shooting method
    fv0 = np.zeros(3,)
    dfv0 = np.array([0., 0., 1.])
    fvHist = np.zeros([Ng, 3])
    def stateRhs(fv, dfv):
        rhs = np.array([fv[1], fv[2], -0.5 * fv[2] * fv[0]])
        J = sps.csr_matrix([[0.0, 1.0, 0.0],
                            [0.0, 0.0, 1.0],
                            [-0.5 * fv[2], 0.0, -0.5 * fv[0]]])
        return rhs, J.dot(dfv)

    Niter = 100
    threshold = 1e-15
    for iter in range(Niter):
        fv = np.copy(fv0)
        dfv = np.copy(dfv0)
        fvHist[0] = np.copy(fv0)
        for n in range(1, Ng):
            k1, dk1 = stateRhs(fv, dfv)
            k2, dk2 = stateRhs(fv + 0.5 * dx * k1, dfv + 0.5 * dx * dk1)
            k3, dk3 = stateRhs(fv + 0.5 * dx * k2, dfv + 0.5 * dx * dk2)
            k4, dk4 = stateRhs(fv + dx * k3, dfv + dx * dk3)
            fv += dx / 6. * (k1 + 2. * k2 + 2. * k3 + k4)
            dfv += dx / 6. * (dk1 + 2. * dk2 + 2. * dk3 + dk4)
            fvHist[n] = np.copy(fv)
        res = fv[1] - 1.0
        if (abs(res) < threshold):
            print(iter, abs(res))
            break
        else:
            fv0[-1] -= res / dfv[1]
    f = fvHist[:, 0]

    temp = np.append(eta[:,np.newaxis].T, fvHist.T, axis=0).T
    # np.savetxt("blasius.txt", temp)
    return temp

def cubicSplineInterp(x, y):
    from scipy.interpolate import splrep, BSpline
    t, c, k = splrep(x, y, k = 3)
    return BSpline(t, c, k)

def initial_condition(g, mach_number=1.0 / 343.0, gamma=1.4, Re = 1.16e6):
    x = g.xyz[0][:,:,:,0]
    y = g.xyz[0][:,:,:,1]

    xmin = np.amin(x)
    delta = np.sqrt((x - xmin) / Re / mach_number + 1.0)
    eta = y / delta
    print(eta[0, 10, -1], eta[-1, 10, -1])
    blasiusTable = solve_blasius(1.5 * np.amax(eta))
    uTable = mach_number * blasiusTable[:, 2]
    vtTable = 0.5 * (blasiusTable[:, 0] * blasiusTable[:, 2] - blasiusTable[:, 1])
    uFunc = cubicSplineInterp(blasiusTable[:, 0], uTable)
    vtFunc = cubicSplineInterp(blasiusTable[:, 0], vtTable)

    u0 = uFunc(eta)
    v0 = vtFunc(eta) / Re / delta
    p = 1.0 / gamma

    s = p3d.Solution().copy_from(g).quiescent(gamma)
    s.q[0][:,:,:,0] = 1.0
    s.q[0][:,:,:,1:4] = 0.0
    s.q[0][:,:,:,1] = u0
    s.q[0][:,:,:,2] = v0
    s.q[0][:,:,:,4] = p
    return s.fromprimitive()

if __name__ == '__main__':
    g = grid([301, 201], mapping_type='geom')
    g.save('DeformingWall2D.xyz')
    s = initial_condition(g, Re = 23175.7)
    s.save('DeformingWall2D.ic.q')
    s.save('DeformingWall2D.target.q')
