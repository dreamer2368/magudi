#!/usr/bin/env python
import sys
import struct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import PLOT3D
from MatplotlibHelper import ApplyTecplotStyle

def drawEdges(ax, grid):
    ax.plot(grid.X[0][:,0,0,0], grid.X[0][:,0,0,1], 'k-')
    ax.plot(grid.X[0][0,:,0,0], grid.X[0][0,:,0,1], 'k-')
    return ax

def plotPressureContours(ax, grid, solutionFile, minLevel = None, 
                         maxLevel = None, nLevels = 31):

    soln = PLOT3D.Solution()
    soln.Import(solutionFile)

    z = (ratioOfSpecificHeats - 1.) * (soln.Q[0][:,:,0,-1] - 
                                       0.5 / soln.Q[0][:,:,0,0] * 
                                       np.sum(soln.Q[0][:,:,0,1:-1] ** 2, axis = -1))

    if minLevel is None: 
        minLevel = z.min()
    if maxLevel is None: 
        maxLevel = z.max()

    ax.contour(grid.X[0][:,:,0,0], grid.X[0][:,:,0,1], z, 
               levels = np.linspace(minLevel, maxLevel, nLevels))

    return ax

if __name__ == '__main__':

    if len(sys.argv) != 2:
        print('Usage: python plotPressureContours.py FILE')
        print('Plots pressure contours from a PLOT3D solution file FILE')
        sys.exit(-1)

    outputPrefix = 'NACA0012'
    ratioOfSpecificHeats = 1.4

    grid = PLOT3D.Grid()
    try:
        grid.Import(outputPrefix + '.xyz')
    except struct.error:
        grid.hasIBLANK = False
        grid.Import(outputPrefix + '.xyz')

    gs = gridspec.GridSpec(1, 1)
    ax = plt.subplot(gs[0,0])
    ApplyTecplotStyle(ax)
    plotPressureContours(ax, grid, sys.argv[1], 
                         minLevel = 1. / ratioOfSpecificHeats - 1e-2,
                         maxLevel = 1. / ratioOfSpecificHeats + 1e-2)
    drawEdges(ax, grid)
    ax.set_xlim([-2., 3.])
    ax.set_ylim([-2., 2.])
    ax.set_aspect('equal')
    plt.show()
