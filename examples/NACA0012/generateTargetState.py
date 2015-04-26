#!/usr/bin/env python
import numpy as np

import PLOT3D

if __name__ == '__main__':

    outputPrefix = 'NACA0012'

    ratioOfSpecificHeats = 1.4
    freeStreamMachNumber = 0.5

    grid = PLOT3D.Grid()
    grid.Import(outputPrefix + '.xyz')

    soln = PLOT3D.Solution()
    soln.CopyFrom(grid)

    soln.Q[0][:,:,0,0] = 1.
    soln.Q[0][:,1:,0,1] = freeStreamMachNumber
    soln.Q[0][:,:,0,-1] = 1. / ratioOfSpecificHeats / (ratioOfSpecificHeats - 1.) + 0.5 / soln.Q[0][:,:,0,0] * soln.Q[0][:,:,0,1] ** 2

    soln.Q[0][-1,:,:,:] = soln.Q[0][0,:,:,:]

    soln.Export(outputPrefix + '.target.q')
