#!/usr/bin/env python
import numpy as np

import PLOT3D

if __name__ == '__main__':

    outputPrefix = 'AcousticMonopole'
    gridSize = [201, 201, 1]

    xMin = -14.
    xMax =  14.
    yMin = -14.
    yMax =  14.

    grid = PLOT3D.Grid()
    grid.nGrids = 1
    grid.SetSize(0, gridSize)
    x = np.linspace(xMin, xMax, gridSize[0])
    y = np.linspace(yMin, yMax, gridSize[1])
    for i in range(x.size):
        grid.X[0][i,:,:,0] = x[i]
    for j in range(y.size):
        grid.X[0][:,j,:,1] = y[j]
    grid.Export(outputPrefix + '.xyz')
