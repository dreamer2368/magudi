#!/usr/bin/env python
import numpy as np

import PLOT3D

def mappingFunction(s, b, c, sigma):
    from scipy.special import erf
    return ((s - 0.5) * (1.0 + 2.0 * b) - b * sigma * (np.exp(- ((s - 0.5 + c) / sigma) ** 2) / np.sqrt(np.pi) + ((s - 0.5 + c) / sigma) * erf((s - 0.5 + c) / sigma) - np.exp(- ((s - 0.5 - c) / sigma) ** 2) / np.sqrt(np.pi) - ((s - 0.5 - c) / sigma) * erf((s - 0.5 - c) / sigma))) / (0.5 + b - b * sigma * (np.exp(- ((0.5 + c) / sigma) ** 2) / np.sqrt(np.pi) + ((0.5 + c) / sigma) * erf((0.5 + c) / sigma) - np.exp(- ((0.5 - c) / sigma) ** 2) / np.sqrt(np.pi) - ((0.5 - c) / sigma) * erf((0.5 - c) / sigma)))

if __name__ == '__main__':

    xMin =  -60.
    xMax =  160.
    yMin = -100.
    yMax =  100.

    gridSize = [960, 640, 1]

    outputPrefix = 'WeiFreundSDML'

    x = mappingFunction(np.linspace(0., 1., gridSize[0]), 80., 0.58, 0.07)
    y = mappingFunction(np.linspace(0., 1., gridSize[1]), 20., 0.62, 0.2)

    grid = PLOT3D.Grid()
    grid.nGrids = 1
    grid.SetSize(0, gridSize)
    for i in range(0, gridSize[0]):
        grid.X[0][i,:,:,0] = xMin + 0.5 * (xMax - xMin) * (1. + x[i])
    for j in range(0, gridSize[1]):
        grid.X[0][:,j,:,1] = yMin + 0.5 * (yMax - yMin) * (1. + y[j])

    grid.Export(outputPrefix + '.xyz')
