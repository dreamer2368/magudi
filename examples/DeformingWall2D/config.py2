#!/usr/bin/env python
import numpy as np
import scipy.sparse as sps
import plot3dnasa as p3d

def mapping_function(x, sigma):
    return np.sinh(sigma * x) / np.sinh(sigma)

def grid(size, Lx, dip_range, mapping_type='sinh'):
    dip_range /= Lx

    from scipy.optimize import fsolve
    x_min =    0. / Lx
    x_max =  180. / Lx
    y_min =    0. / Lx
    y_max =   30. / Lx
    z_min =  -40. / Lx
    z_max =   40. / Lx
    dy_min = 0.016
    num_uniform = 7

    x = np.linspace(x_min, x_max, size[0] + 1)[:-1]
    if mapping_type == 'sinh':
        sigma = fsolve(lambda x: (
                (y_max - y_min - num_uniform * dy_min) * mapping_function(
                    1. / (size[1] - 1.), x) - dy_min) ** 2, 2.)
        y = np.append(np.linspace(y_min, y_min + dy_min * num_uniform,
                                  num_uniform + 1), y_min + dy_min *
                      num_uniform + (y_max - y_min - num_uniform * dy_min) *
                      mapping_function(np.linspace(0., 1., size[1] -
                                                   num_uniform), sigma)[1:])
    else:
        sigma = fsolve(lambda x: (y_max - y_min) / dy_min - num_uniform + 1 -
                       (x ** (size[1] - num_uniform) - 1.) /
                       (x - 1.), 1.02)
        print 100. * (sigma - 1.)
        y = np.append([0.], np.cumsum(
                [dy_min if r < num_uniform - 1
                 else dy_min * sigma ** (r - num_uniform + 1)
                 for r in range(size[1] - 1)]))

    # leftIdx, rightIdx = p3d.find_extents(x, dip_range[0], dip_range[1])
    dummy, depthIdx = p3d.find_extents(y, 0.0, dip_range[2])
    # assert(leftIdx >= 1)
    # assert(rightIdx <= size[0])
    assert(depthIdx <= size[1])
    size[1] += depthIdx - 1
    y = np.append(-y[(depthIdx-1):0:-1], y)

    # sizes = [size, [rightIdx - leftIdx + 1, depthIdx]]
    # print(sizes)
    g = p3d.Grid().set_size(size, True)
    # z = np.linspace(z_min, z_max, g.size[0,2] + 1)[:-1]
    for i in range(x.size):
        g.xyz[0][i,:,:,0] = x[i]
    for j in range(y.size):
        g.xyz[0][:,j,:,1] = y[j]
    # for k in range(z.size):
    #     g.xyz[0][:,:,k,2] = z[k]

    # for i in range(rightIdx - leftIdx + 1):
    #     g.xyz[1][i,:,:,0] = x[i + leftIdx - 1]
    # for j in range(depthIdx):
    #     g.xyz[1][:,-1-j,:,1] = -y[j]

    patchRange = int(1.5 * depthIdx)
    formatStr = '  {:<20} {:<21} {:>4d} {:>7d}' + 6 * ' {:>4d}'
    print(formatStr.format('immersed1', 'IMMERSED_BOUNDARY',
                           1, 0, 1, -1, 1, patchRange, 1, -1))
    # print(formatStr.format('interface.S1', 'SAT_BLOCK_INTERFACE',
    #                        1, 2, leftIdx, rightIdx, 1, 1, 1, -1))
    # print(formatStr.format('interface.N2', 'SAT_BLOCK_INTERFACE',
    #                        2, -2, 1, -1, -1, -1, 1, -1))
    # print(formatStr.format('immersed1', 'IMMERSED_BOUNDARY',
    #                        1, 0, leftIdx, rightIdx, 1, depthIdx, 1, -1))
    # print(formatStr.format('immersed2', 'IMMERSED_BOUNDARY',
    #                        2, 0, 1, -1, 1, -1, 1, -1))
    # print(formatStr.format('flatPlate.S1', 'SAT_ISOTHERMAL_WALL',
    #                        1, 2, 1, leftIdx - 1, 1, 1, 1, -1))
    # print(formatStr.format('flatPlate.N1', 'SAT_ISOTHERMAL_WALL',
    #                        1, 2, rightIdx + 1, -1, 1, 1, 1, -1))
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

    mask = (eta < 0.0)
    u0[mask] = 0.0
    v0[mask] = 0.0

    s = p3d.Solution().copy_from(g).quiescent(gamma)
    s.q[0][:,:,:,0] = 1.0
    s.q[0][:,:,:,1:4] = 0.0
    s.q[0][:,:,:,1] = u0
    s.q[0][:,:,:,2] = v0
    s.q[0][:,:,:,4] = p
    return s.fromprimitive()

def initial_condition_from_exp(g, expFilename, Lx = 5.0, a0 = 343.0, gamma = 1.4, Re = 1.1588e+05):
    x = g.xyz[0][:,:,:,0]
    y = g.xyz[0][:,:,:,1]

    expData = np.loadtxt(expFilename)
    ymax = expData[-1, 0] / Lx
    uFunc = cubicSplineInterp(expData[:, 0], expData[:, 1])
    vtFunc = cubicSplineInterp(expData[:, 0], expData[:, 2])

    u0 = uFunc(y * Lx) / a0
    v0 = vtFunc(y * Lx) / a0
    p = 1.0 / gamma

    mask = (y < 0.0)
    u0[mask] = 0.0
    v0[mask] = 0.0

    u0[u0 <= 0.0] = 0.0
    v0[v0 <= 0.0] = 0.0
    u0[y > ymax] = expData[-1, 1] / a0
    v0[y > ymax] = expData[-1, 2] / a0

    s = p3d.Solution().copy_from(g).quiescent(gamma)
    s.q[0][:,:,:,0] = 1.0
    s.q[0][:,:,:,1:4] = 0.0
    s.q[0][:,:,:,1] = u0
    s.q[0][:,:,:,2] = v0
    s.q[0][:,:,:,4] = p
    return s.fromprimitive()

if __name__ == '__main__':
    Lx = 5.0        # mm
    a0 = 343.0      # m/s
    gamma = 1.4
    nu0 = 1.48e-5   # m2/s
    Re = Lx * 1e-3 * a0 / nu0
    print("Reynolds number: %.5E" % Re)
    expFile = 'UV.4.5v_0msC.txt'

    dip_range = np.array([0.0, 15.0, 7.5])
    g = grid([501, 201], Lx, dip_range, mapping_type='geom')
    g.save('DeformingWall2D.xyz')
    s = initial_condition_from_exp(g, expFile, Lx, a0, gamma, Re)
    # s = initial_condition(g, mach_number = 0.1, Re = 23175.7)
    s.save('DeformingWall2D.ic.q')
    s.save('DeformingWall2D.target.q')
