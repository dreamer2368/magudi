#ini!/usr/bin/env python
import numpy as np
import plot3dnasa as p3d

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

def mapping_function(s, b, c, sigma):
    from scipy.special import erf
    return ((s - 0.5) * (1.0 + 2.0 * b) - b * sigma *
            (np.exp(- ((s - 0.5 + c) / sigma) ** 2) / np.sqrt(np.pi) +
             ((s - 0.5 + c) / sigma) * erf((s - 0.5 + c) / sigma) -
             np.exp(- ((s - 0.5 - c) / sigma) ** 2) / np.sqrt(np.pi) -
             ((s - 0.5 - c) / sigma) * erf((s - 0.5 - c) / sigma))) / \
        (0.5 + b - b * sigma * (np.exp(- ((0.5 + c) / sigma) ** 2) /
                                np.sqrt(np.pi) + ((0.5 + c) / sigma) *
                                erf((0.5 + c) / sigma) -
                                np.exp(- ((0.5 - c) / sigma) ** 2) /
                                np.sqrt(np.pi) - ((0.5 - c) / sigma) *
                                erf((0.5 - c) / sigma)))

def grid(size,xMin,xMax,yMin,yMax):
		g = p3d.Grid().set_size(size, True)
		x = np.linspace(xMin,xMax,size[0])
		s = np.linspace(0., 1., size[1])
		background_y=np.linspace(yMin,yMax,size[1])
		xWallStart = -5.
		xWallEnd   = 5.
		wallHeight = wallProfile((x - xWallStart) / (xWallEnd - xWallStart), 0.05 , 2e-4, 10)
	
		for i in range(size[0]):
			g.xyz[0][i,:,0,0] = x[i]

		for j in range(size[1]):
				yo=background_y[j]
				width=1./4.
			
				if s[j] <= width:
					xi=1.*((width-(width-s[j]))/width)**0.5
				else:
					xi=1.

				#print "j={j},y={y},xi={xi}".format(j=j,y=yo,xi=xi)	
				for k in range(size[0]):
					allege_y=s[j] * yMax + (1. - s[j]) * (yMin + wallHeight[k])
					g.xyz[0][k,j,0,1] = allege_y - xi*(allege_y-yo)

		return g

def target_state(g, u1, u2, S, gamma=1.4):
    s = p3d.Solution().copy_from(g).quiescent(gamma)
    x = g.xyz[0][:,:,:,0]
    s.q[0][:,:,:,1] = u2 + 0.5 * (u1 - u2) * (
        1. + np.tanh(2. * g.xyz[0][:,:,:,1] / \
                     (1. + S * np.where(x > 0., x, np.zeros_like(x)))))
    return s.fromprimitive(gamma)

def initial_condition(g, u1, u2, gamma=1.4):
    s = p3d.Solution().copy_from(g).quiescent(gamma)
    s.q[0][:,:,:,1] = u2 + 0.5 * (u1 - u2) * (
        1. + np.tanh(2. * g.xyz[0][:,:,:,1]))
    return s.fromprimitive(gamma)

def target_mollifier(g,x_min,x_max,y_min,y_max):
    f = p3d.Function().copy_from(g)
    f.f[0].fill(1.)
    n = f.get_size(0)
    for j in range(n[1]):
        f.f[0][:,j,0,0] *= p3d.tanh_support(
            g.xyz[0][:,j,0,0], x_min, x_max, 40., 0.2)    
    for i in range(n[0]):
        f.f[0][i,:,0,0] *= p3d.cubic_bspline_support(
            g.xyz[0][i,:,0,1], y_min, y_max)
    return f

def control_mollifier(g,x_min,x_max,y_min,y_max):
    f = p3d.Function().copy_from(g)
    f.f[0].fill(1.)
    n = f.get_size(0)
    for j in range(n[1]):
        f.f[0][:,j,0,0] *= p3d.cubic_bspline_support(
            g.xyz[0][:,j,0,0], x_min, x_max)
    for i in range(n[0]):
        f.f[0][i,:,0,0] *= p3d.cubic_bspline_support(
            g.xyz[0][i,:,0,1], y_min, y_max)
    return f

def mean_pressure(s):
    f = p3d.Function().copy_from(s)
    f.f[0][:,:,:,0] = s.toprimitive().q[0][:,:,:,4]
    return f

if __name__ == '__main__':

		outputPrefix = 'WavyWallAcousticMonopole'
        
		xMin = -10.
		xMax =  10.
		yMin = -10.
		yMax =  10.

		g = grid([201,201],xMin,xMax,yMin,yMax)
		g.save(outputPrefix+'.xyz')
		initial_condition(g, u1=0.0, u2=0.0).save(outputPrefix+'.ic.q')
		target_state(g, u1=0.0, u2=0.0, S=0.05).save(outputPrefix+'.target.q')

		xMinTarget = 5.
		xMaxTarget = 7.
		yMinTarget = -9.
		yMaxTarget = 8.
		target_mollifier(g,xMinTarget,xMaxTarget,yMinTarget,yMaxTarget).save(outputPrefix+'.target_mollifier.f')

		xMinControl =  -6.
		xMaxControl =  -4.
		yMinControl = -8.5
		yMaxControl =  -6.5
		control_mollifier(g,xMinControl,xMaxControl,yMinControl,yMaxControl).save(outputPrefix+'.control_mollifier.f')

		#IMPLEMENTING THE BOUNDARY CONDITION FILE

		#SPONGE INDICES
		left_sponge=np.argmin(abs(g.xyz[0][:,0,0,0] - (xMin+2))) + 1
		right_sponge=np.argmin(abs(g.xyz[0][:,0,0,0] -(xMax-2)) ) + 1
		top_sponge=np.argmin(abs(g.xyz[0][0,:,0,1] - (yMax-2))) + 1

		#TARGET INDICES
		left_target=np.argmin(abs(g.xyz[0][:,0,0,0] -(xMinTarget)) ) + 1
		right_target=np.argmin(abs(g.xyz[0][:,0,0,0] - (xMaxTarget))) + 1
		bottom_target=np.argmin(abs(g.xyz[0][0,:,0,1] - (yMinTarget))) + 1
		top_target=np.argmin(abs(g.xyz[0][0,:,0,1] - (yMaxTarget))) + 1

		#CONTROL INDICES
		left_c=np.argmin(abs(g.xyz[0][:,0,0,0] - (xMinControl))) + 1
		right_c=np.argmin(abs(g.xyz[0][:,0,0,0] - (xMaxControl))) + 1
		bottom_c=np.argmin(abs(g.xyz[0][0,:,0,1] - (yMinControl))) + 1
		top_c=np.argmin(abs(g.xyz[0][0,:,0,1] - (yMaxControl))) + 1

		with open(outputPrefix+'_bc.dat', "wt") as fp:
				fp.write("# Name                 Type                  Grid normDir iMin iMax jMin jMax kMin kMax\n")
				fp.write("# ==================== ===================== ==== ======= ==== ==== ==== ==== ==== ====\n")
				fp.write("inflow               SAT_FAR_FIELD            1       1    1    1    1   -1    1   -1\n")
				fp.write("inflowSponge         SPONGE                   1       1    1  %i    1   -1    1   -1\n"% (left_sponge))
				fp.write("outflow              SAT_FAR_FIELD            1      -1   -1   -1    1   -1    1   -1\n")
				fp.write("outflowSponge        SPONGE                   1      -1   %i   -1    1   -1    1   -1\n"% (right_sponge))
				fp.write("wall.S               SAT_ISOTHERMAL_WALL            1       2    1   -1    1    1    1   -1\n")
				fp.write("farfield.N           SAT_FAR_FIELD            1      -2    1   -1   -1   -1    1   -1\n")
				fp.write("sponge.N             SPONGE                   1      -2    1   -1  %i   -1    1   -1\n"% (top_sponge))
				fp.write("targetRegion         COST_TARGET              1       0  %i %i  %i  %i    1   -1\n"%(left_target,right_target,bottom_target,top_target))
				fp.write("controlRegion        ACTUATOR                 1       0  %i  %i %i %i 1   -1\n"%(left_c,right_c,bottom_c,top_c))
