#!/usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as plt

import PLOT3D

def plot(ax, dataFile, *args, **kwargs):

    a = np.loadtxt(dataFile, usecols = [1, 2])
    ax.plot(a[:,0], a[:,1], *args, **kwargs.pop('plot_kwargs', dict()))

    return ax

if __name__ == '__main__':

    from textwrap import fill

    if len(sys.argv) != 2:
        print(fill('Usage: python plotInstantaneousCost.py FILE', 80) + '\n')
        print(fill('Plots the instantaneous cost functional from a text file FILE', 80) + '\n')
        print(fill('FILE is assumed to have three columns with the second and third column representing sampling times and the values of the instantaneous cost functional, respectively.', 80))
        sys.exit(-1)

    outputPrefix = 'WavyWall'
    
    ax = plt.subplot(111)
    plot(ax, sys.argv[1], 'k-', plot_kwargs = dict())
    plt.show()
