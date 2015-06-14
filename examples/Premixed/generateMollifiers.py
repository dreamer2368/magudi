#!/usr/bin/env python
import numpy as np

import PLOT3D

def cubic_bspline_support(x, x_min, x_max):
    from scipy.signal import cubic
    imin = np.argmin(np.abs(x - x_min))
    imax = np.argmin(np.abs(x - x_max))
    assert imax - imin + 1 >= 3
    return cubic(4. * (x - x[imin]) / (x[imax] - x[imin]) - 2.)

def tanh_support(x, x_min, x_max, sigma, xi):
    imin = np.argmin(np.abs(x - x_min))
    imax = np.argmin(np.abs(x - x_max))
    assert imax > imin
    f = lambda x: \
        np.tanh(sigma * (x + 1. - 0.5 * xi)) - \
        np.tanh(sigma * (x - 1. + 0.5 * xi))
    y = f(2. * (x - x[imin]) / (x[imax] - x[imin]) - 1.)
    return y - y.min()

def find_extents(x, x_min, x_max):
    return np.where(x >= x_min)[0][0], np.where(x <= x_max)[0][-1] + 2

if __name__ == '__main__':

    outputPrefix = 'Premixed'

    grid = PLOT3D.Grid()
    grid.Import(outputPrefix + '.xyz')

    xMinTarget = -1.
    xMaxTarget =  1.
    yMinTarget = -10.
    yMaxTarget =  10.

    iMinTarget, iMaxTarget = find_extents(grid.X[0][:,0,0,0], xMinTarget, xMaxTarget)
    jMinTarget, jMaxTarget = find_extents(grid.X[0][0,:,0,1], yMinTarget, yMaxTarget)
    print ('  {:<20} {:<21} {:>4d} {:>7d}' + 6 * ' {:>4d}').format('targetRegion', 'COST_TARGET', 1, 0, iMinTarget, iMaxTarget, jMinTarget, jMaxTarget, 1, -1)

    xMinControl =  1.
    xMaxControl =  5.
    yMinControl = -2.
    yMaxControl =  2.

    iMinControl, iMaxControl = find_extents(grid.X[0][:,0,0,0], xMinControl, xMaxControl)
    jMinControl, jMaxControl = find_extents(grid.X[0][0,:,0,1], yMinControl, yMaxControl)
    print ('  {:<20} {:<21} {:>4d} {:>7d}' + 6 * ' {:>4d}').format('controlRegion', 'ACTUATOR', 1, 0, iMinControl, iMaxControl, jMinControl, jMaxControl, 1, -1)

    func = PLOT3D.Function()
    func.nComponents = 1
    func.CopyFrom(grid)

    for l in range(grid.nGrids):
        gridSize = grid.GetSize(l)
        func.F[l][:,:,:,0] = 1.
        for i in range(gridSize[0]):
            func.F[l][i,:,:,0] *= tanh_support(grid.X[l][i,:,:,1], yMinTarget, yMaxTarget, 40., 0.2)
        for j in range(gridSize[1]):
            func.F[l][:,j,:,0] *= cubic_bspline_support(grid.X[l][:,j,:,0], xMinTarget, xMaxTarget)
    func.Export(outputPrefix + '.target_mollifier.f')

    for l in range(grid.nGrids):
        gridSize = grid.GetSize(l)
        func.F[l][:,:,:,0] = 1.
        for i in range(gridSize[0]):
            func.F[l][i,:,:,0] *= cubic_bspline_support(grid.X[l][i,:,:,1], yMinControl, yMaxControl)
        for j in range(gridSize[1]):
            func.F[l][:,j,:,0] *= cubic_bspline_support(grid.X[l][:,j,:,0], xMinControl, xMaxControl)
    func.Export(outputPrefix + '.control_mollifier.f')
