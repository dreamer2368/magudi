#!/usr/bin/env python
import numpy as np

import PLOT3D

def CubicBsplineSupport(x, xMin, xMax):

    from scipy.signal import cubic

    iMin = np.argmin(abs(x - xMin))
    iMax = np.argmin(abs(x - xMax))
    assert iMax - iMin + 1 >= 3

    return cubic(4. * (x - x[iMin]) / (x[iMax] - x[iMin]) - 2.)

def HyperbolicTangentSupport(x, xMin, xMax, steepness, rampFraction):

    iMin = np.argmin(abs(x - xMin))
    iMax = np.argmin(abs(x - xMax))
    assert iMax > iMin

    f = lambda x: \
        np.tanh(steepness * (x + 1. - 0.5 * rampFraction)) - \
        np.tanh(steepness * (x - 1. + 0.5 * rampFraction))

    y = f(2. * (x - x[iMin]) / (x[iMax] - x[iMin]) - 1.)

    return y - y.min()

def findExtents(x, xMin, xMax):

    # Fortran indexing
    iMin = 1
    iMax = x.size

    while x[iMin-1] < xMin:
        iMin += 1
    iMin -= 1

    while x[iMax-1] > xMax:
        iMax -= 1
    iMax += 1

    assert iMin >= 1 and iMin <= x.size
    assert iMax >= 1 and iMax <= x.size

    return iMin, iMax

if __name__ == '__main__':

    outputPrefix = 'ReactiveJet'

    grid = PLOT3D.Grid()
    grid.Import(outputPrefix + '.xyz')

    xMinTarget =  20.
    xMaxTarget =  36.
    yMinTarget = -8.
    yMaxTarget =  8.

    iMinTarget, iMaxTarget = findExtents(grid.X[0][:,0,0,0], xMinTarget, xMaxTarget)
    jMinTarget, jMaxTarget = findExtents(grid.X[0][0,:,0,1], yMinTarget, yMaxTarget)
    print ('  {:<20} {:<21} {:>4d} {:>7d}' + 6 * ' {:>4d}').format('targetRegion', 'COST_TARGET', 1, 0, iMinTarget, iMaxTarget, jMinTarget, jMaxTarget, 1, -1)

    xMinControl =  4.
    xMaxControl =  12.
    yMinControl = -4.
    yMaxControl =  4.

    iMinControl, iMaxControl = findExtents(grid.X[0][:,0,0,0], xMinControl, xMaxControl)
    jMinControl, jMaxControl = findExtents(grid.X[0][0,:,0,1], yMinControl, yMaxControl)
    print ('  {:<20} {:<21} {:>4d} {:>7d}' + 6 * ' {:>4d}').format('controlRegion', 'ACTUATOR', 1, 0, iMinControl, iMaxControl, jMinControl, jMaxControl, 1, -1)

    func = PLOT3D.Function()
    func.nComponents = 1
    func.CopyFrom(grid)

    for l in range(grid.nGrids):
        gridSize = grid.GetSize(l)
        func.F[l][:,:,:,0] = 1.
        for i in range(gridSize[0]):
            func.F[l][i,:,:,0] *= CubicBsplineSupport(grid.X[l][i,:,:,1], yMinTarget, yMaxTarget)
        for j in range(gridSize[1]):
            func.F[l][:,j,:,0] *= CubicBsplineSupport(grid.X[l][:,j,:,0], xMinTarget, xMaxTarget)
    func.Export(outputPrefix + '.target_mollifier.f')

    for l in range(grid.nGrids):
        gridSize = grid.GetSize(l)
        func.F[l][:,:,:,0] = 1.
        for i in range(gridSize[0]):
            func.F[l][i,:,:,0] *= CubicBsplineSupport(grid.X[l][i,:,:,1], yMinControl, yMaxControl)
        for j in range(gridSize[1]):
            func.F[l][:,j,:,0] *= CubicBsplineSupport(grid.X[l][:,j,:,0], xMinControl, xMaxControl)
    func.Export(outputPrefix + '.control_mollifier.f')
