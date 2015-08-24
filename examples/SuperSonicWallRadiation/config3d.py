#ini!/usr/bin/env python
import numpy as np
import plot3dnasa as p3d

def ShapeMollifierFunction(xNormalized):
	MollifierFunction = lambda x: 0.5*(np.tanh(40. * (x - 0.1)) - np.tanh(40. * (x -0.6))) #np.where(np.logical_or(x < 0., x > 60), 0., np.tanh(40. * (x - 10)) - np.tanh(40. * (x - 50)))
	Mollifier = MollifierFunction(xNormalized)
	#Mollifier /= 0.8 #np.trapz(Mollifier, xNormalized) # normalized to have unit area

	return Mollifier

def wallProfile(size,xNormalized,zNormalized, amplitude, nModes):


	from numpy.random import rand
	from numpy.random import seed
	from numpy import meshgrid

	wallHeight = np.zeros((len(xNormalized), len(zNormalized)))
	#shapeMollifier = ShapeMollifierFunction(xNormalized)

	seed(1985)
	controlParameters = rand(nModes)
	phaseAngles = 2. * np.pi * rand(nModes)
	xwaveNumbers = 2. * np.pi * np.arange(1, nModes+1)
	zwaveNumbers=np.zeros((2*nModes+1))
	for i in range(len(zwaveNumbers)):
		zwaveNumbers[i]=2.*np.pi*(-nModes+i)

	phases=np.zeros((nModes,nModes*2+1))
	for m in range(len(zwaveNumbers)):
		for n in range(nModes):
			phases[n,m]=2.*np.pi*rand()

	for i in range(len(xNormalized)):
		for j in range(len(zNormalized)):
			mollx=0.5*(np.tanh(40. * (xNormalized[i] - 0.1)) - np.tanh(40. * (xNormalized[i] -0.6)))
			mollz=0.5*(np.tanh(40. * (zNormalized[j] - 0.05)) - np.tanh(40. * (zNormalized[j] -0.95)))
			for m in range(len(zwaveNumbers)):
				for n in range(nModes):
					wallHeight[i,j] += amplitude *mollx*mollz* np.cos(zwaveNumbers[m] * zNormalized[j] + xwaveNumbers[n]*xNormalized[i]+phases[n,m])

	wallHeight *=1.	
	return wallHeight

def grid(size,xMin,xMax,yMin,yMax,zMin,zMax):
		g = p3d.Grid().set_size(size, True)
		x = np.linspace(xMin,xMax,size[0])
		z = np.linspace(zMin,zMax,size[2])
		s = np.linspace(0., 1., size[1])
		background_y=np.linspace(yMin,yMax,size[1])
		#0.0006324555325
		wallHeight = wallProfile(size,np.linspace(0., 1., size[0]),np.linspace(0.,1.,size[2]),0.001,6)

		for i in range(size[2]):
			g.xyz[0][:,:,i,2] = z[i]	
		for i in range(size[0]):
			g.xyz[0][i,:,:,0] = x[i]

		for j in range(size[1]):
				yo=background_y[j]
				width=1./5.
			
				if s[j] <= width:
					xi=np.tanh(5.*s[j]/width)
				else:
					xi=1.

				for l in range(size[2]):
					for k in range(size[0]):
						allege_y=s[j] * yMax + (1. - s[j]) * (yMin + wallHeight[k,l])
						g.xyz[0][k,j,l,1] = allege_y - xi*(allege_y-yo)
		return g


def ambient_state(g, gamma=1.4):
     s = p3d.Solution().copy_from(g).quiescent(gamma)
     return s.fromprimitive(gamma)

def target_state(g, u1, u2, S,yCenter, gamma=1.4):
    s = p3d.Solution().copy_from(g).quiescent(gamma)
    x = g.xyz[0][:,:,:,0]
    s.q[0][:,:,:,1] = u2 + 0.5 * (u1 - u2) * (
        1. + np.tanh(2. * (g.xyz[0][:,:,:,1]-yCenter) / \
                     (1. + S * np.where((x) > 0.,x, np.zeros_like(x)))))
    return s.fromprimitive(gamma)

def initial_condition(g, u1, u2,S,yCenter, gamma=1.4):
	s = p3d.Solution().copy_from(g).quiescent(gamma)
	s.q[0][:,:,:,1] = u2 + 0.5 * (u1 - u2) * (
	   1. + np.tanh(2. * (g.xyz[0][:,:,:,1]-yCenter)))
	return s.fromprimitive(gamma)

def target_mollifier(g,x_min,x_max,y_min,y_max):
    f = p3d.Function().copy_from(g)
    f.f[0].fill(1.)
    n = f.get_size(0)
    for j in range(n[1]):
        f.f[0][:,j,0,0] *= p3d.tanh_support(
            g.xyz[0][:,j,0,0], x_min, x_max, 40., 0.2)    
    for i in range(n[0]):
		f.f[0][i,:,0,0] *= p3d.tanh_support(
        	g.xyz[0][i,:,0,1], y_min, y_max,40., 0.2)	
    return f

def control_mollifier(g,x_min,x_max,y_min,y_max):
	n = g.get_size(0)
	x = np.linspace(x_min,x_max,n[0])
	shapeMollifier=ShapeMollifierFunction(np.linspace(0.,1.,n[0]))

	f = p3d.Function().copy_from(g)
	f.f[0].fill(0.)

	for i in range(n[0]):
	   f.f[0][i,0,0,0]=shapeMollifier[i]
	return f

def mean_pressure(s):
    f = p3d.Function().copy_from(s)
    f.f[0][:,:,:,0] = s.toprimitive().q[0][:,:,:,4]
    return f

def ambient_pressure(s):
     f = p3d.Function().copy_from(s)
     f.f[0][:,:,:,0] = s.toprimitive().q[0][:,:,:,4]
     return f


if __name__ == '__main__':

		outputPrefix = 'SuperSonicWallRadiation'
        
		xMin = 0.
		xMax = 2.
		yMin = 0.
		yMax = 1.
		zMin=0.
		zMax=1.5

		g = grid([201,101,151],xMin,xMax,yMin,yMax,zMin,zMax)
		g.save(outputPrefix+'.xyz')

		yShearMax=-10.
		initial_condition(g, u1=1.25, u2=0.0,S=0.0,yCenter=yShearMax).save(outputPrefix+'.ic.q')
		target_state(g, u1=1.25, u2=0.0, S=0.0,yCenter=yShearMax).save(outputPrefix+'.target.q')

		ambient_state(g,gamma=1.4).save(outputPrefix+'.ambient.q')
		ambient_pressure(p3d.fromfile(outputPrefix+'.ambient.q')).save(outputPrefix+'.ambient_pressure.f')

		xMinTarget = 0.3
		xMaxTarget = 0.7
		yMinTarget = 0.4
		yMaxTarget = 0.5
		target_mollifier(g,xMinTarget,xMaxTarget,yMinTarget,yMaxTarget).save(outputPrefix+'.target_mollifier.f')

		xMinControl =  xMin
		xMaxControl =  xMax
		yMinControl =  yMin
		yMaxControl =  yMin
		control_mollifier(g,xMin,xMax,yMin,yMax).save(outputPrefix+'.control_mollifier.f')

		#IMPLEMENTING THE BOUNDARY CONDITION FILE

		 #SPONGE INDICES
		right_sponge=np.argmin(abs(g.xyz[0][:,0,0,0] -(xMax-0.1)) ) + 1
		left_sponge=np.argmin(abs(g.xyz[0][:,0,0,0] -(xMin+0.1)) ) + 1
		top_sponge=np.argmin(abs(g.xyz[0][0,:,0,1] - (yMax-0.1))) + 1
		front_sponge=np.argmin(abs(g.xyz[0][0,0,:,2] - (zMin+0.1))) + 1
		back_sponge=np.argmin(abs(g.xyz[0][0,0,:,2] - (zMax-0.1))) + 1

		#TARGET INDICES
		left_target=np.argmin(abs(g.xyz[0][:,0,0,0] -(xMinTarget)) ) + 1
		right_target=np.argmin(abs(g.xyz[0][:,0,0,0] - (xMaxTarget))) + 1
		bottom_target=np.argmin(abs(g.xyz[0][0,:,0,1] - (yMinTarget))) + 1
		top_target=np.argmin(abs(g.xyz[0][0,:,0,1] - (yMaxTarget))) + 1

		#CONTROL INDICES
		left_c=np.argmin(abs(g.xyz[0][:,0,0,0] - (xMinControl))) + 1
		right_c=np.argmin(abs(g.xyz[0][:,0,0,0] - (xMaxControl))) + 1
		bottom_c=np.argmin(abs(g.xyz[0][0,:,0,1] - (yMinControl))) + 1
		top_c=np.argmin(abs(g.xyz[0][0,:,0,1] - (yMinControl))) + 1


		with open(outputPrefix+'_bc.dat', "wt") as fp:
				fp.write("# Name                 Type                  Grid normDir iMin iMax jMin jMax kMin kMax\n")
				fp.write("# ==================== ===================== ==== ======= ==== ==== ==== ==== ==== ====\n")
				fp.write("inflow               SAT_FAR_FIELD            1       1    1    1    1   -1    1   -1\n")
				fp.write("inflowSponge         SPONGE            1       1    1    %i    1   -1    1   -1\n"% (left_sponge))
				fp.write("outflow              SAT_FAR_FIELD            1      -1   -1   -1    1   -1    1   -1\n")
				fp.write("outflowSponge        SPONGE                   1      -1   %i   -1    1   -1    1   -1\n"% (right_sponge))
				fp.write("outflowFront              SAT_FAR_FIELD            1      3   1   -1    1   -1    1   1\n")
				fp.write("outflowFrontSponge        SPONGE                   1      3   1   -1    1   -1    1   %i\n"% (front_sponge))
				fp.write("outflowBack              SAT_FAR_FIELD            1      -3   1   -1    1   -1    -1   -1\n")
				fp.write("outflowBackSponge        SPONGE                   1      -3   1   -1    1   -1    %i   -1\n"% (back_sponge))
				fp.write("wall.S               SAT_SLIP_WALL            1       2    1   -1    1    1    1   -1\n")
				fp.write("farfield.N           SAT_FAR_FIELD            1      -2    1   -1   -1   -1    1   -1\n")
				fp.write("sponge.N             SPONGE                   1      -2    1   -1  %i   -1    1   -1\n"% (top_sponge))
    				#fp.write("targetRegion         COST_TARGET              1       0  %i %i  %i  %i    1   -1\n"%(left_target,right_target,bottom_target,top_target))
				#fp.write("controlRegion.p1        ACTUATOR                 1       2  %i  %i %i %i 1   -1\n"%(left_c,right_c,bottom_c,top_c))
