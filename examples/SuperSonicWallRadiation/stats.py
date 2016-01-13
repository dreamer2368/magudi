#!/usr/bin/env python
import numpy as np
import plot3dnasa as p3d


def ambient_state(g, gamma=1.4):
	s = p3d.Solution().copy_from(g).quiescent(gamma)
	return s.fromprimitive(gamma)

def target_state(g, u1, u2, S,yCenter, gamma=1.4):
	s = p3d.Solution().copy_from(g).quiescent(gamma)
	x = g.xyz[0][:,:,:,0]
	s.q[0][:,:,:,1] = u2 + 0.5 * (u1 - u2) * (1. + np.tanh(2. * (g.xyz[0][:,:,:,1]-yCenter) / \
		(1. + S * np.where((x) > 0.,x, np.zeros_like(x)))))
	return s.fromprimitive(gamma)

def initial_condition(g, u1, u2,S,yCenter, gamma=1.4):
	s = p3d.Solution().copy_from(g).quiescent(gamma)
	s.q[0][:,:,:,1] = u2 + 0.5 * (u1 - u2) * (1. + np.tanh(2. * (g.xyz[0][:,:,:,1]-yCenter)))
	return s.fromprimitive(gamma)

def ambient_pressure(s):
	f = p3d.Function().copy_from(s)
	f.f[0][:,:,:,0] = s.toprimitive().q[0][:,:,:,4]
	return f

def read_sol(external_grid,external_solution,gamma=1.4):
	ge = p3d.Grid(external_grid)
	xe = ge.load().xyz[0][:,0,0,0]
	ye = ge.load().xyz[0][0,:,0,1]
	ze = ge.load().xyz[0][0,0,:,2]
	se = p3d.fromfile(external_solution)
	p = p3d.Solution().copy_from(ge).quiescent(gamma)
	p=se.toprimitive()

	f = open('output.dat', 'w')
	f.write('y,yrms,p_sk,p_rms/pavg,v_rms\n')

	for j in range(0,50): 
		#RMS OF QUANTITIES
		p_rms=0.
		v_rms=0.
		rho_rms=0.
		y_rms=0.
		p_sk=0.
		rho_sk=0.
	
		pavg=0.
		vavg=0.
		rhoavg=0.
		yavg=0.

		atJ=j
		M=1.5
		beta=(M**2.-1)**0.5
		imax=xe.size
		IntegralX=7.46

		yProbe=ge.xyz[0][0,atJ,0,1]
		xmin=0.+yProbe*beta
		xmax=ge.xyz[0][imax-1,0,0,0]

		imin=np.argmin(abs(ge.xyz[0][:,0,0,0] - (IntegralX+yProbe*beta))) + 1
		imax=np.argmin(abs(ge.xyz[0][:,0,0,0] - (IntegralX*2+yProbe*beta))) + 1	

		N=(imax-imin+1)*ze.size

		for k in range(ze.size):
			for i in range(imin,imax+1):
				pavg=pavg+p.q[0][i,atJ,k,4]
				vavg=vavg+p.q[0][i,atJ,k,2]
				rhoavg=rhoavg+p.q[0][i,atJ,k,0]
				yavg=yavg+ge.xyz[0][i,atJ,k,1]

		pavg=pavg/N
		vavg=vavg/N
		rhoavg=rhoavg/N
		yavg=yavg/N

		for k in range(ze.size):
			for i in range(imin,imax+1):
				p_rms=p_rms+(p.q[0][i,atJ,k,4]-pavg)**2.
				v_rms=v_rms+(p.q[0][i,atJ,k,2]-vavg)**2.
				rho_rms=rho_rms+(p.q[0][i,atJ,k,0]-rhoavg)**2.
				y_rms=y_rms+(ge.xyz[0][i,atJ,k,1]-yavg)**2.
				p_sk=p_sk+(p.q[0][i,atJ,k,4]-pavg)**3.
				rho_sk=rho_sk+(p.q[0][i,atJ,k,0]-rhoavg)**3.

		p_rms=(p_rms/N)**0.5
		v_rms=(v_rms/N)**0.5
		rho_rms=(rho_rms/N)**0.5
		y_rms=(y_rms/N)**0.5
		p_sk=p_sk/N
		rho_sk=rho_sk/N

		p_sk=p_sk/(p_rms)**3
		rho_sk=rho_sk/(rho_rms)**3
		f.write("%f,%f,%f,%f,%f\n"%(yProbe,y_rms**0.5,p_sk,p_rms/pavg,v_rms))	

		print yProbe,y_rms**0.5,p_sk,p_rms/pavg,v_rms,0.239,0.283,0.325


	return 

if __name__ == '__main__':
		read_sol(external_grid='SuperSonicWallRadiation.xyz',external_solution='SuperSonicWallRadiation-00003000.q')
