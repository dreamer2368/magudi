#!/usr/bin/env python
import numpy as np
import plot3dnasa as p3d

def grid(size,alpha=1.0):
    x_min = 0.0
    x_max = alpha
    y_min = 0.0
    y_max = 1.0
    g = p3d.Grid().set_size(size, True)
    x = np.linspace(x_min, x_max, g.size[0,0] + 1)[:-1]
    y = np.linspace(y_min, y_max, g.size[0,1] + 1)[:-1]
    g.xyz[0][:,:,0,:2] = np.transpose(np.meshgrid(x, y))
    return g

def initial_condition(g, mach_number=0.1, gamma=1.4):
    y = g.xyz[0][:,:,:,1]
    x = g.xyz[0][:,:,:,0]

    temp = np.zeros_like(g.xyz[0][:,:,:,0])
    temp = 1.0 / (gamma-1)

    s = p3d.Solution().copy_from(g).quiescent(gamma)
    s.q[0][:,:,:,0] = 1.0
    s.q[0][:,:,:,1:4] = 0.0
    s.q[0][:,:,:,2] = mach_number * np.cos(2.0*np.pi*(x-0.1*np.sin(2.0*np.pi*y)))
    s.q[0][:,:,:,4] = 1.0 / gamma
    return s.fromprimitive()

def random_solution(g,time=0.0,timestep=0):

    gamma = 1.4
    s = p3d.Solution().copy_from(g).quiescent(gamma)
    n = s.get_size(0)
    s.q[0][:,:,0,0] = ( 2.0 * np.random.rand(n[0],n[1]) - 1.0 ) * 1.0e-2 + 1.0
    s.q[0][:,:,0,1] = ( 2.0 * np.random.rand(n[0],n[1]) - 1.0 ) * 1.0e-2 + 0.0
    s.q[0][:,:,0,2] = ( 2.0 * np.random.rand(n[0],n[1]) - 1.0 ) * 1.0e-2 + 0.1
    s.q[0][:,:,0,4] = ( 2.0 * np.random.rand(n[0],n[1]) - 1.0 ) * 1.0e-2 + 1./(gamma-1.)
    s.q[0][:,:,0,4] *= (gamma-1.)/gamma * s.q[0][:,:,0,0]

    s.time = time
    s._format.aux_header[0] = timestep

    return s.fromprimitive(gamma)

if __name__ == '__main__':
    g = grid([256,256],alpha=1.0)
    g.save('KolmogorovFlow.xyz')
    initial_condition(g,mach_number=0.2).save('KolmogorovFlow.ic.q')
    #random_solution(g).save('ChannelFlow.ic.q')
