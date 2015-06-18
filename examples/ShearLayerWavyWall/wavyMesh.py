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

def Wei_mappingFunction(s, b, c, sigma):
    from scipy.special import erf
    return ((s - 0.5) * (1.0 + 2.0 * b) - b * sigma * (np.exp(- ((s - 0.5 + c) / sigma) ** 2) / np.sqrt(np.pi) + ((s - 0.5 + c) / sigma) * erf((s - 0.5 + c) / sigma) - np.exp(- ((s - 0.5 - c) / sigma) ** 2) / np.sqrt(np.pi) - ((s - 0.5 - c) / sigma) * erf((s - 0.5 - c) / sigma))) / (0.5 + b - b * sigma * (np.exp(- ((0.5 + c) / sigma) ** 2) / np.sqrt(np.pi) + ((0.5 + c) / sigma) * erf((0.5 + c) / sigma) - np.exp(- ((0.5 - c) / sigma) ** 2) / np.sqrt(np.pi) - ((0.5 - c) / sigma) * erf((0.5 - c) / sigma)))

if __name__ == '__main__':

	outputPrefix = 'WeiFreundSDML'
	xMin =  -60.
	xMax =  160.
	yMin = -20.
	yMax =  100.

	xMinTarget =   0.
	xMaxTarget = 100.
	yMinTarget = 65.
	yMaxTarget = 75.

	xMinControl =  1.
	xMaxControl =  7.
	yMinControl = -3.
	yMaxControl =  3.

	xMinExcite=-15.
	xMaxExcite=-5.
	yMinExcite=-2.
	yMaxExcite=2.

	gridSize = [400,400, 1]

	grid = PLOT3D.Grid()
	grid.scalarType = np.float64
	grid.nGrids = 1
	grid.SetSize(0, gridSize)

	gx = Wei_mappingFunction(np.linspace(0., 1., gridSize[0]), 80., 0.58, 0.07) #as done in Wei Thesis pg 14
	x = xMin + 0.5 * (xMax - xMin) * (1. + gx)

	s = np.linspace(0., 1., gridSize[1])

	xWallStart = 5.
	xWallEnd   = 60.
	wallHeight = wallProfile((x - xWallStart) / (xWallEnd - xWallStart), 0.25, 2e-4, 8)

	for i in range(gridSize[0]):
		grid.X[0][i,:,0,0] = x[i]
	for j in range(gridSize[1]):
		for k in range(gridSize[0]):
			grid.X[0][k,j,0,1] = s[j] * yMax + (1. - s[j]) * (yMin + wallHeight[k])
	grid.Export(outputPrefix + '.xyz')

	#SPONGE INDICES
	left_sponge=np.argmin(abs(grid.X[0][:,0,0,0] + (xMin+60))) + 1
	right_sponge=np.argmin(abs(grid.X[0][:,0,0,0] -(xMax-60)) ) + 1
	print "left_sponge= ",left_sponge,grid.X[0][left_sponge,0,0,0]
	print "right_sponge= ",right_sponge,grid.X[0][right_sponge,0,0,0]

	top_sponge=np.argmin(abs(grid.X[0][0,:,0,1] - (yMax-20))) + 1
	print "top_sponge= ",top_sponge,grid.X[0][0,top_sponge,0,1]

	#TARGET INDICES
	right_target=right_sponge
	left_target=left_sponge
	bottom_target=np.argmin(abs(grid.X[0][0,:,0,1] - (yMinTarget))) + 1
	top_target=np.argmin(abs(grid.X[0][0,:,0,1] - (yMaxTarget))) + 1
	print "bot target= ",bottom_target
	print "top target= ",top_target

	#CONTROL INDICES
	left_c=np.argmin(abs(grid.X[0][:,0,0,0] - (xMinControl))) + 1
	right_c=np.argmin(abs(grid.X[0][:,0,0,0] - (xMaxControl))) + 1
	bottom_c=np.argmin(abs(grid.X[0][0,:,0,1] - (yMinControl))) + 1
	top_c=np.argmin(abs(grid.X[0][0,:,0,1] - (yMaxControl))) + 1

	#EXCITATION INDICES
	left_ex=np.argmin(abs(grid.X[0][:,0,0,0] - (xMinExcite))) + 1
	right_ex=np.argmin(abs(grid.X[0][:,0,0,0] - (xMaxExcite))) + 1
	bottom_ex=np.argmin(abs(grid.X[0][0,:,0,1] - (yMinExcite))) + 1
	top_ex=np.argmin(abs(grid.X[0][0,:,0,1] - (yMaxExcite))) + 1

	with open(outputPrefix+'_bc.dat', "wt") as fp: 
		fp.write("# Name                 Type                  Grid normDir iMin iMax jMin jMax kMin kMax\n")
		fp.write("# ==================== ===================== ==== ======= ==== ==== ==== ==== ==== ====\n")
		fp.write("inflow               SAT_FAR_FIELD            1       1    1    1    1   -1    1   -1\n") 
		fp.write("inflowSponge         SPONGE                   1       1    1  %i    1   -1    1   -1\n"% (left_sponge)) 
		fp.write("outflow              SAT_FAR_FIELD            1      -1   -1   -1    1   -1    1   -1\n") 
		fp.write("outflowSponge        SPONGE                   1      -1   %i   -1    1   -1    1   -1\n"% (right_sponge))
		fp.write("wall.S               SAT_SLIP_WALL            1       2    1   -1    1    1    1   -1\n")  
		fp.write("farfield.N           SAT_FAR_FIELD            1      -2    1   -1   -1   -1    1   -1\n") 
		fp.write("sponge.N             SPONGE                   1      -2    1   -1  %i   -1    1   -1\n"% (top_sponge))
		fp.write("excitationSupport    SOLENOIDAL_EXCITATION    1       0   %i  %i  %i   %i    1   -1\n"%(left_ex,right_ex,bottom_ex,top_ex))
		fp.write("targetRegion         COST_TARGET              1       0  %i %i  %i  %i    1   -1\n"%(left_target,right_target,bottom_target,top_target))
       		fp.write("controlRegion        ACTUATOR                 1       0  %i  %i %i %i 1   -1\n"%(left_c,right_c,bottom_c,top_c))
