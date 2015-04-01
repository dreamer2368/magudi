#!/usr/bin/env python
import os
import sys
sys.path.append(os.getcwd())

import numpy as np
from scipy.interpolate import griddata
from scipy.ndimage.filters import gaussian_filter1d, uniform_filter1d

import PLOT3D

if __name__ == '__main__':

    outputPrefix = 'NACA0012'
    ratioOfSpecificHeats = 1.4
    freeStreamDensity = 1.225
    freeStreamPressure = 101325.
    gasConstant = 8.314462175 * 1000. / 28.966
    specificHeatAtConstantPressure = 1006.43
    specificHeatAtConstantVolume = specificHeatAtConstantPressure - gasConstant
    ratioOfSpecificHeats = specificHeatAtConstantPressure / specificHeatAtConstantVolume
    freeStreamSpeedOfSound = np.sqrt(ratioOfSpecificHeats * freeStreamPressure / freeStreamDensity)
    freeStreamMachNumber = 0.2

    a = np.loadtxt(outputPrefix + '_steady_ascii.txt', skiprows = 1, usecols = range(1, 7))
    x = a[:,0]
    y = a[:,1]
    density = a[:,2] / freeStreamDensity
    xVelocity = a[:,3] / freeStreamSpeedOfSound
    yVelocity = a[:,4] / freeStreamSpeedOfSound
    pressure = (a[:,5] + freeStreamPressure) / freeStreamPressure / ratioOfSpecificHeats

    grid = PLOT3D.Grid()
    grid.hasIBLANK = False
    grid.Import(outputPrefix + '_PLOT3D.x')
    X = np.copy(grid.X[0][:,:,0,0])
    Y = np.copy(grid.X[0][:,:,0,1])
    grid.Update(0)
    grid.X[0][:,:,0,0] = X
    grid.X[0][:,:,0,1] = Y
    grid.Export(outputPrefix + '.xyz')

    soln = PLOT3D.Solution()
    soln.CopyFrom(grid)    

    # Linear interpolation from unstructured grid used by Fluent.
    soln.Q[0][:,:,0,0] = 1.
    soln.Q[0][:,:,0,1] = griddata((x, y), xVelocity, (X, Y), method = 'linear', fill_value = freeStreamMachNumber)
    soln.Q[0][:,:,0,2] = griddata((x, y), yVelocity, (X, Y), method = 'linear', fill_value = 0.)
    
    # Filter velocity field.

    nFilterRepeat = 50

    showProgress = True
    try:
        from progressbar import ProgressBar
    except:
        showProgress = False

    if showProgress is True:
        progressBar = ProgressBar(maxval = nFilterRepeat)
        progressBar.start()

    for i in range(nFilterRepeat):
        soln.Q[0][:-1,:,0,1:3] = uniform_filter1d(soln.Q[0][:-1,:,0,1:3], size = 3, axis = 0, mode = 'wrap')
        soln.Q[0][:-1,:,0,1:3] = uniform_filter1d(soln.Q[0][:-1,:,0,1:3], size = 3, axis = 1, mode = 'nearest')
        if showProgress is True:
            progressBar.update(i)
    
    if showProgress is True:
        progressBar.finish()

    pressure = 1. / ratioOfSpecificHeats + 0.5 * freeStreamMachNumber ** 2 * (1. - np.sum(soln.Q[0][:,:,0,1:3] ** 2, axis = -1) / freeStreamMachNumber ** 2)
    soln.Q[0][:,:,0,-1] = pressure / (ratioOfSpecificHeats - 1.) + 0.5 / soln.Q[0][:,:,0,0] * np.sum(soln.Q[0][:,:,0,1:3] ** 2, axis = -1)

    # Adjust all variables on overlapping plane.
    soln.Q[0][-1,:,:,:] = soln.Q[0][0,:,:,:]

    soln.Export(outputPrefix + '.ic.q')

    # Target state.
    soln.Q[0][:,:,0,0] = 1.
    soln.Q[0][:,:,0,1] = freeStreamMachNumber
    soln.Q[0][:,:,0,-1] = 1. / ratioOfSpecificHeats / (ratioOfSpecificHeats - 1.) + 0.5 / soln.Q[0][:,:,0,0] * soln.Q[0][:,:,0,1] ** 2
    soln.Export(outputPrefix + '.target.q')
