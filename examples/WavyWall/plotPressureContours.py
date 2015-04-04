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

    ratioOfSpecificHeats = kwargs.pop('ratioOfSpecificHeats', 1.4)

    z = (ratioOfSpecificHeats - 1.) * (soln.Q[0][:,:,0,-1] - 0.5 / soln.Q[0][:,:,0,0] * np.sum(soln.Q[0][:,:,0,1:-1] ** 2, axis = -1))

    ax.contourf(grid.X[0][:,:,0,0], grid.X[0][:,:,0,1], z, *args, **kwargs.pop('plot_kwargs', dict()))

    return ax

if __name__ == '__main__':

    from textwrap import fill

    if len(sys.argv) != 2:
        print(fill('Usage: python plotPressureContours.py FILE', 80) + '\n')
        print(fill('Plots pressure contours for a two-dimensional flow from a PLOT3D solution file FILE', 80) + '\n')
        print(fill('If FILE is a multi-block PLOT3D file, only the data corresponding to the first block is plotted.', 80))
        sys.exit(-1)

    outputPrefix = 'WavyWall'

    ratioOfSpecificHeats = 1.4

    grid = PLOT3D.Grid()
    try:
        grid.Import(outputPrefix + '.xyz')
    except struct.error:
        grid.hasIBLANK = False
        grid.Import(outputPrefix + '.xyz')

    ax = plt.subplot(111)
    plot(ax, grid, sys.argv[1], ratioOfSpecificHeats = ratioOfSpecificHeats, plot_kwargs = dict(extend = 'both', levels = np.linspace(1. / ratioOfSpecificHeats - 1e-3, 1. / ratioOfSpecificHeats + 1e-3, 31)))
    drawEdges(ax, grid)
    ax.set_xlim([-10., 10.])
    ax.set_ylim([-11., 10.])
    ax.set_aspect('equal')
    plt.show()
