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

if __name__ == '__main__':

    outputPrefix = 'AcousticMonopole'

    xMinTarget =  -1.
    xMaxTarget =   1.
    yMinTarget = -10.
    yMaxTarget =  10.

    xMinControl =  1.
    xMaxControl =  5.
    yMinControl = -2.
    yMaxControl =  2.

    grid = PLOT3D.Grid()
    grid.scalarType = np.float64
    grid.Import(outputPrefix + '.xyz')

    func = PLOT3D.Function()
    func.nComponents = 1
    func.CopyFrom(grid)

    for l in range(grid.nGrids):
        gridSize = grid.GetSize(l)
        assert gridSize[2] == 1
        func.F[l][:,:,0,0] = 1.
        for i in range(gridSize[0]):
            func.F[l][i,:,0,0] *= HyperbolicTangentSupport(grid.X[l][i,:,0,1], yMinTarget, yMaxTarget, 40., 0.2)
        for j in range(gridSize[1]):
            func.F[l][:,j,0,0] *= CubicBsplineSupport(grid.X[l][:,j,0,0], xMinTarget, xMaxTarget)
    func.Export(outputPrefix + '.target_mollifier.f')

    for l in range(grid.nGrids):
        gridSize = grid.GetSize(l)
        assert gridSize[2] == 1
        func.F[l][:,:,0,0] = 1.
        for i in range(gridSize[0]):
            func.F[l][i,:,0,0] *= CubicBsplineSupport(grid.X[l][i,:,0,1], yMinControl, yMaxControl)
        for j in range(gridSize[1]):
            func.F[l][:,j,0,0] *= CubicBsplineSupport(grid.X[l][:,j,0,0], xMinControl, xMaxControl)
    func.Export(outputPrefix + '.control_mollifier.f')
