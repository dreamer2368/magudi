#!/usr/bin/env python
import numpy as np
from pylab import *
from scipy.optimize import fsolve, minimize
import scipy.fftpack

import PLOT3D

def airfoilProfile(x, chordLength, thicknessRatio):
    return 5. * thicknessRatio * chordLength * (0.2969 * np.sqrt(x / chordLength) - 0.1260 * (x / chordLength) - 0.3516 * (x / chordLength) ** 2 + 0.2843 * (x / chordLength) ** 3 - 0.1015 * (x / chordLength) ** 4)

def airfoilProfileSlope(x, chordLength, thicknessRatio):
    return 5. * thicknessRatio * (0.14845 / np.sqrt(x / chordLength) - 0.1260 - 0.7032 * (x / chordLength) + 0.8529 * (x / chordLength) ** 2 - 0.406 * (x / chordLength) ** 3)

def KarmanTrefftzTransform(z, zSingularPoint1, zSingularPoint2, zetaSingularPoint1, zetaSingularPoint2, includedAngle):
    transformationExponent = np.pi / (2. * np.pi - includedAngle)
    return zetaSingularPoint2 + (zetaSingularPoint2 - zetaSingularPoint1) / (((z - zSingularPoint1) / (z - zSingularPoint2)) ** transformationExponent - 1.)

def KarmanTrefftzInverseTransform(zeta, zSingularPoint1, zSingularPoint2, zetaSingularPoint1, zetaSingularPoint2, includedAngle):
    transformationExponent = 1. / (np.pi / (2. * np.pi - includedAngle))
    return zSingularPoint2 + (zSingularPoint2 - zSingularPoint1) / (((zeta - zetaSingularPoint1) / (zeta - zetaSingularPoint2)) ** transformationExponent - 1.)

def KarmanTrefftzInverseTransformDerivative(zeta, zSingularPoint1, zSingularPoint2, zetaSingularPoint1, zetaSingularPoint2, includedAngle):
    transformationExponent = 1. / (np.pi / (2. * np.pi - includedAngle))
    return - transformationExponent * (zSingularPoint2 - zSingularPoint1) / (((zeta - zetaSingularPoint1) / (zeta - zetaSingularPoint2)) ** transformationExponent - 1.) ** 2 * ((zeta - zetaSingularPoint1) / (zeta - zetaSingularPoint2)) ** (transformationExponent - 1.) * (zetaSingularPoint1 - zetaSingularPoint2) / (zeta - zetaSingularPoint2) ** 2

if __name__ == '__main__':

    outputPrefix = 'NACA0012'
    gridSize = [129, 129, 1]

    chordLength = 10.
    thicknessRatio = 0.12

    ratioOfSpecificHeats = 1.4
    freeStreamMachNumber = 0.2
    angleOfAttack = - 2. * np.pi / 180.

    mappingTolerance = 1e-14

    xTrailingEdge, = fsolve(airfoilProfile, chordLength, args = (chordLength, thicknessRatio))
    leadingEdgeCurvature = 0.5 / chordLength * (5. * chordLength * thicknessRatio * 0.2969) ** 2
    includedAngle = 2. * np.arctan(abs(airfoilProfileSlope(xTrailingEdge, chordLength, thicknessRatio)))    

    # Physical coordinate: z
    zSingularPoint1 = xTrailingEdge + 0.j
    zSingularPoint2 = 0.5 * leadingEdgeCurvature + 0.j

    # Mapped coordinate: zeta    
    zetaSingularPoint1 = 0.77043505 + 0.j
    zetaSingularPoint2 = 0.24642903 + 0.j
    zetaLeadingEdge = KarmanTrefftzTransform(0. + 0.j, zSingularPoint1, zSingularPoint2, zetaSingularPoint1, zetaSingularPoint2, includedAngle)
    zetaCenter = 0.5 * (zetaLeadingEdge + zetaSingularPoint1)
    
    assert (gridSize[0] - 1) % 2 == 0
    assert gridSize[0] >= gridSize[1]

    theta = np.linspace(0., np.pi, (gridSize[0] + 1) / 2)
    radiusInMappedPlane = np.empty_like(theta)
    radiusInMappedPlane.fill(np.abs(zetaCenter - zetaSingularPoint1))

    for i, theta_i in enumerate(theta[1:-1]):
        for j in range(100):
            zetaAirfoil = zetaCenter + radiusInMappedPlane[i] * np.exp(1.j * theta_i)
            zAirfoil = KarmanTrefftzInverseTransform(zetaAirfoil, zSingularPoint1, zSingularPoint2, zetaSingularPoint1, zetaSingularPoint2, includedAngle)
            yExact = airfoilProfile(zAirfoil.real, chordLength, thicknessRatio)
            if np.abs(zAirfoil.imag - yExact) < mappingTolerance:
                break
            zAirfoil = zAirfoil.real + 1j * yExact
            zetaAirfoil = KarmanTrefftzTransform(zAirfoil, zSingularPoint1, zSingularPoint2, zetaSingularPoint1, zetaSingularPoint2, includedAngle)
            radiusInMappedPlane[i] = np.abs(zetaAirfoil - zetaCenter)

    theta = np.linspace(0., 2. * np.pi, gridSize[0])
    radiusInMappedPlane = np.append(radiusInMappedPlane, radiusInMappedPlane[-2::-1])
    zetaAirfoil = zetaCenter + radiusInMappedPlane * np.exp(1.j * theta)

    concentricCirclesRadii = np.exp(2. * np.pi * np.linspace(0., 1., gridSize[0]))
    quasiCircleArcLength = np.sum(np.sqrt(scipy.fftpack.diff(radiusInMappedPlane[:-1]) ** 2 + radiusInMappedPlane[:-1] ** 2) * np.diff(theta))
    concentricCirclesRadii *= quasiCircleArcLength / (2. * np.pi)

    zeta = np.empty([gridSize[0], gridSize[0]], dtype = zetaAirfoil.dtype)
    for i in range(gridSize[0]):
        zeta[:,i] = zetaCenter + (radiusInMappedPlane * (concentricCirclesRadii[-1] - concentricCirclesRadii[0]) + concentricCirclesRadii[-1] * (concentricCirclesRadii[i] - concentricCirclesRadii[0])) / (concentricCirclesRadii[-1] - concentricCirclesRadii[0]) * np.exp(1.j * theta)

    grid = PLOT3D.Grid()
    grid.nGrids = 1
    grid.SetSize(0, gridSize)

    z = KarmanTrefftzInverseTransform(zeta, zSingularPoint1, zSingularPoint2, zetaSingularPoint1, zetaSingularPoint2, includedAngle)
    for i in range(gridSize[0]):
        x = 0.5 * (z[:(gridSize[0]+1)/2,i].real + z[-1:(gridSize[0]-1)/2-1:-1,i].real)
        y = 0.5 * (z[:(gridSize[0]+1)/2,i].imag - z[-1:(gridSize[0]-1)/2-1:-1,i].imag)
        if i == 0:
            x[-1] = 0.
        y[-1] = 0.
        x = np.append(x,  x[-2::-1])
        y = np.append(y, -y[-2::-1])
        if i < gridSize[1]:
            grid.X[0][:,i,0,0] = x
            grid.X[0][:,i,0,1] = y
        z[:,i] = x + 1.j * y

    # Make the grid right-handed.
    grid.X[0][:,:,:,1] = - grid.X[0][:,:,:,1]
    grid.X[0][-1,:,:,:] = grid.X[0][0,:,:,:]
        
    grid.Export(outputPrefix + '.xyz')

    uMappedPlane = np.empty_like(zeta)
    for i in range(gridSize[0]):
        circulationGamma = 4. * np.pi * freeStreamMachNumber * radiusInMappedPlane[0] * np.sin(angleOfAttack)
        uMappedPlane[i,:] = freeStreamMachNumber * (np.exp(-1.j * angleOfAttack) - radiusInMappedPlane[0] ** 2 / (zeta[i,:] - zetaCenter) ** 2 * np.exp(1.j * angleOfAttack)) + 1.j * circulationGamma / (2. * np.pi * (zeta[i,:] - zetaCenter))

    error = uMappedPlane[:,0] * (zeta[:,0] - zetaCenter) / radiusInMappedPlane[0]
    print np.abs(error.real).max()
    uPhysicalPlane = uMappedPlane / KarmanTrefftzInverseTransformDerivative(zeta, zSingularPoint1, zSingularPoint2, zetaSingularPoint1, zetaSingularPoint2, includedAngle)

    soln = PLOT3D.Solution()
    soln.CopyFrom(grid)

    soln.Q[0][:,:,0,0] = 1.
    for i in range(gridSize[1]):
        soln.Q[0][:,i,0,1] = uPhysicalPlane[:,i].real
        soln.Q[0][:,i,0,2] = uPhysicalPlane[:,i].imag

    # Hack for scaling correctly
    soln.Q[0][:,:,0,1:3] = freeStreamMachNumber * soln.Q[0][:,:,0,1:3] / np.mean(soln.Q[0][:,-1,0,1])

    pressure = 1. / ratioOfSpecificHeats + 0.5 * (freeStreamMachNumber ** 2 - (soln.Q[0][:,:,0,1] ** 2 + soln.Q[0][:,:,0,2] ** 2))

    soln.Q[0][:,:,0,-1] = pressure / (ratioOfSpecificHeats - 1.) + 0.5 / soln.Q[0][:,:,0,0] * (soln.Q[0][:,:,0,1] ** 2 + soln.Q[0][:,:,0,2] ** 2)

    soln.Export(outputPrefix + '.target.q')

