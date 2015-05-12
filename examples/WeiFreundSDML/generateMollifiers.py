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

    outputPrefix = 'WeiFreundSDML'

    xMinTarget =   0.
    xMaxTarget = 100.
    yMinTarget = -75.
    yMaxTarget = -65.

    xMinControl =  1.
    xMaxControl =  7.
    yMinControl = -3.
    yMaxControl =  3.

    grid = PLOT3D.Grid()
    grid.Import(outputPrefix + '.xyz')

    func = PLOT3D.Function()
    func.nComponents = 1
    func.CopyFrom(grid)

    for l in range(grid.nGrids):
        gridSize = grid.GetSize(l)
        func.F[l][:,:,:,0] = 1.
        for i in range(gridSize[0]):
            func.F[l][i,:,:,0] *= CubicBsplineSupport(grid.X[l][i,:,:,1], yMinTarget, yMaxTarget)
        for j in range(gridSize[1]):
            func.F[l][:,j,:,0] *= HyperbolicTangentSupport(grid.X[l][:,j,:,0], xMinTarget, xMaxTarget, 80., 0.2)
    func.Export(outputPrefix + '.target_mollifier.f')

    for l in range(grid.nGrids):
        gridSize = grid.GetSize(l)
        func.F[l][:,:,:,0] = 1.
        for i in range(gridSize[0]):
            func.F[l][i,:,:,0] *= CubicBsplineSupport(grid.X[l][i,:,:,1], yMinControl, yMaxControl)
        for j in range(gridSize[1]):
            func.F[l][:,j,:,0] *= CubicBsplineSupport(grid.X[l][:,j,:,0], xMinControl, xMaxControl)
    func.Export(outputPrefix + '.control_mollifier.f')
