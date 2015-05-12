#!/usr/bin/env python
import numpy as np

import PLOT3D

if __name__ == '__main__':

    outputPrefix = 'WeiFreundSDML'

    grid = PLOT3D.Grid()
    grid.Import(outputPrefix + '.xyz')

    x = grid.X[0][:,:,:,0]
    y = grid.X[0][:,:,:,1]

    lowerFluidVelocity = 0.2
    upperFluidVelocity = 0.9

    ratioOfSpecificHeats = 1.4

    soln = PLOT3D.Solution()
    soln.CopyFrom(grid)

    soln.Q[0][:,:,:,0] = 1.
    soln.Q[0][:,:,:,1] = lowerFluidVelocity + 0.5 * (upperFluidVelocity - lowerFluidVelocity) * (1. + np.tanh(2. * y))
    soln.Q[0][:,:,:,-1] = 1. / ratioOfSpecificHeats / (ratioOfSpecificHeats - 1.) + 0.5 / soln.Q[0][:,:,:,0] * np.sum(soln.Q[0][:,:,:,1:-1] ** 2, axis = -1)

    soln.Export(outputPrefix + '.ic.q')
