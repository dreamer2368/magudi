#!/usr/bin/env python
import sys
import struct
import numpy as np
import matplotlib.pyplot as plt

import PLOT3D

def drawEdges(ax, grid):
    ax.plot(grid.X[0][:,0,0,0], grid.X[0][:,0,0,1], 'k-')
    ax.plot(grid.X[0][0,:,0,0], grid.X[0][0,:,0,1], 'k-')
    return ax

def plot(ax, grid, solutionFile, *args, **kwargs):

    soln = PLOT3D.Solution()
    soln.CopyFrom(grid)
    soln.Import(solutionFile)

    z = soln.Q[0][:,:,0,-1]

    ax.contourf(grid.X[0][:,:,0,0], grid.X[0][:,:,0,1], z, **kwargs.pop('plot_kwargs', dict()))

    return ax

if __name__ == '__main__':

    from textwrap import fill

    if len(sys.argv) != 2:
        print(fill('Usage: python plotAdjointPressureContours.py FILE', 80) + '\n')
        print(fill('Plots adjoint pressure contours for a two-dimensional flow from a PLOT3D solution file FILE', 80))
        print(fill('If FILE is a multi-block PLOT3D file, only the data corresponding to the first block is plotted.', 80))
        sys.exit(-1)

    outputPrefix = 'WavyWall'

    grid = PLOT3D.Grid()
    try:
        grid.Import(outputPrefix + '.xyz')
    except struct.error:
        grid.hasIBLANK = False
        grid.Import(outputPrefix + '.xyz')

    ax = plt.subplot(111)
    plot(ax, grid, sys.argv[1], plot_kwargs = dict(extend = 'both', levels = np.linspace(-1e-5, 1e-5, 30)))
    drawEdges(ax, grid)
    ax.set_xlim([-10., 10.])
    ax.set_ylim([-11., 10.])
    ax.set_aspect('equal')
    plt.show()
