#!/usr/bin/env python
import numpy as np

import PLOT3D

def wallProfile(xNormalized, amplitude, gaussianFactor, nModes):

    from numpy.random import rand

    wallHeight = np.zeros_like(xNormalized)
    
    shapeMollifierFunction = lambda x: np.where(np.logical_or(x < 0., x > 1.), 0., np.tanh(40. * (x - 0.1)) - np.tanh(40. * (x - 0.9)))
    shapeMollifier = shapeMollifierFunction(xNormalized)
    shapeMollifier /= np.trapz(shapeMollifier, xNormalized) # normalized to have unit area

    controlParameters = rand(nModes)
    phaseAngles = 2. * np.pi * rand(nModes)
    waveNumbers = 2. * np.pi * np.arange(1, nModes + 1)
    
    for i in range(nModes):
        wallHeight += amplitude * shapeMollifier * controlParameters[i] * np.exp(- gaussianFactor * waveNumbers[i] ** 2) * np.cos(waveNumbers[i] * xNormalized + phaseAngles[i])

    return wallHeight

if __name__ == '__main__':

    outputPrefix = 'WavyWall'
    gridSize = [201, 201, 1]

    xMin = -14.
    xMax =  14.
    yMin = -10.
    yMax =  14.

    grid = PLOT3D.Grid()
    grid.scalarType = np.float64
    grid.nGrids = 1
    grid.SetSize(0, gridSize)
    x = np.linspace(xMin, xMax, gridSize[0])
    s = np.linspace(0., 1., gridSize[1])

    xWallStart = -10.
    xWallEnd   = +10.
    wallHeight = wallProfile((x - xWallStart) / (xWallEnd - xWallStart), 0.1, 2e-4, 6)

    for i in range(gridSize[0]):
        grid.X[0][i,:,0,0] = x[i]
    for j in range(gridSize[1]):
        grid.X[0][:,j,0,1] = s[j] * yMax + (1. - s[j]) * (yMin + wallHeight)
    grid.Export(outputPrefix + '.xyz')
