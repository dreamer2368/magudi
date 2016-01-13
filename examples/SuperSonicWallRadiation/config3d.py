#ini!/usr/bin/env python
import numpy as np
import plot3dnasa as p3d

def ShapeMollifierFunction(xNormalized):
	MollifierFunction = lambda x: 0.5*(np.tanh(25. * (x - 0.05)) - np.tanh(25. * (x -1.5))) #np.where(np.logical_or(x < 0., x > 60), 0., np.tanh(40. * (x - 10)) - np.tanh(40. * (x - 50)))
	Mollifier = MollifierFunction(xNormalized)
	#Mollifier /= 0.8 #np.trapz(Mollifier, xNormalized) # normalized to have unit area

	return Mollifier


def wallProfile(size,xNormalized,zNormalized,amplitude,kx,kz,xMin,xMax):

	wallHeight = np.zeros((len(xNormalized), len(zNormalized)))

	for i in range(len(xNormalized)):
		mollx=1.
		#mollx=0.5*(np.tanh(1. * (xNormalized[i] - 0.25*2.*np.pi/kx)) - np.tanh(1. * (xNormalized[i] -(xMax+0.25*2*np.pi/kx))))
		#mollx=0.5*(np.tanh(1. * (xNormalized[i] - (-10.))) - np.tanh(1. * (xNormalized[i] -(xMax-np.pi/kx))))
		for j in range(len(zNormalized)):
			wallHeight[i,j] +=  mollx*amplitude*(np.sin(kz*zNormalized[j] + kx*xNormalized[i]) +np.sin(-kz*zNormalized[j] + kx*xNormalized[i]))
	
	return wallHeight

def grid(size,xMin,xMax,yMin,yMax,zMin,zMax,kx,kz,epsilon):

		g = p3d.Grid().set_size(size, True)
		x = np.linspace(xMin,xMax,size[0])
		z = np.linspace(zMin,zMax,size[2])
		s = np.linspace(0., 1., size[1])
		background_y=np.linspace(yMin,yMax,size[1])
		
		wallHeight = wallProfile(size,x,z,epsilon,kx,kz,xMin,xMax)

		for i in range(z.size):
			g.xyz[0][:,:,i,2] = z[i]	
		for i in range(size[0]):
			g.xyz[0][i,:,:,0] = x[i]

		for j in range(size[1]):
				yo=background_y[j]
				width=1./2.
			
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
   
		Nx=301
		Nz=151
		Ny=151
 
		Integral_x=9.51E+00
		Integral_z=2.61E+00
		numXModes=4
		numZModes=int(numXModes*Integral_x/Nx*Nz/Integral_z)
	
		Lpx=numXModes*Integral_x
		Lpz=numZModes*Integral_z   

		xMin = 0.
		xMax = Lpx/Nx * (Nx-1)
		yMin = 0.
		yMax = xMax*0.5
		zMin=0.
		zMax=Lpz/Nz * (Nz-1)

		kx=2.*np.pi/Integral_x
		kz=2.*np.pi/Integral_z
		epsilon=0.1

		right_sponge_width=0.05*xMax

		print 'xMax=',xMax
		print 'zMax=',zMax
		print 'Lpz=',zMax+zMax/(Nz-1)
		print 'Lpx=',xMax+xMax/(Nx-1)
		print 'sin(kz*zMax)=',np.sin(kz*zMax)
		print 'sin(kz*Lpz)=',np.sin(kz*Lpz)
		print 'sin(kz*(zMax+zMax/(Nz-1)))',np.sin(kz*(zMax+zMax/(Nz-1)))
		print 'dx,dy,dz',xMax/(Nx-1),yMax/(Ny-1),zMax/(Nz-1)
		#output linspace to file to ensure we are periodic

		g = grid([Nx,Ny,Nz],xMin,xMax,yMin,yMax,zMin,zMax,kx,kz,epsilon)
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
		#target_mollifier(g,xMinTarget,xMaxTarget,yMinTarget,yMaxTarget).save(outputPrefix+'.target_mollifier.f')

		xMinControl =  xMin
		xMaxControl =  xMax
		yMinControl =  yMin
		yMaxControl =  yMin
		#control_mollifier(g,xMin,xMax,yMin,yMax).save(outputPrefix+'.control_mollifier.f')

		#IMPLEMENTING THE BOUNDARY CONDITION FILE

		 #SPONGE INDICES
		right_sponge_width=0.25*2*np.pi/kx
		right_sponge=np.argmin(abs(g.xyz[0][:,0,0,0] -(xMax-right_sponge_width)) ) + 1
		top_sponge=np.argmin(abs(g.xyz[0][0,:,0,1] - (yMax-right_sponge_width))) + 1

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
				#fp.write("inflow               SAT_FAR_FIELD            1       1    1    1    1   -1    1   -1\n")
				fp.write("outflow              SAT_FAR_FIELD            1      -1   -1   -1    1   -1    1   -1\n")
				fp.write("wall.S               SAT_SLIP_WALL            1       2    1   -1    1    1    1   -1\n")
				fp.write("farfield.N           SAT_FAR_FIELD            1      -2    1   -1   -1   -1    1   -1\n")
				#fp.write("outflowSponge        SPONGE                   1      -1   %i   -1    1   -1    1   -1\n"% (right_sponge))
				fp.write("sponge.N             SPONGE                   1      -2    1   -1  %i   -1    1   -1\n"% (top_sponge))
