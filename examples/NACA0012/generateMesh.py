#!/usr/bin/env python
import numpy as np

import PLOT3D

if __name__ == '__main__':

    outputPrefix = 'NACA0012'

    grid = PLOT3D.Grid()
    grid.hasIBLANK = False
    grid.Import('Gridgen_' + outputPrefix + '_PLOT3D.x')
    gridSize = grid.GetSize(0)
    X = np.copy(grid.X[0][:,:,0,0])
    Y = np.copy(grid.X[0][:,:,0,1])
    print gridSize

    grid = PLOT3D.Grid()
    grid.nGrids = 1
    grid.SetSize(0, gridSize)
    grid.X[0][:,:,0,0] = X
    grid.X[0][:,:,0,1] = Y
    grid.Export(outputPrefix + '.xyz')
