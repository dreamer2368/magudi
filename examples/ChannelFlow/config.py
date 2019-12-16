#!/usr/bin/env python
import numpy as np
import plot3dnasa as p3d

def mapping_function(x, sigma):
    return np.sinh(sigma * x) / np.sinh(sigma)

def grid(size, eta=0.97):
    from scipy.optimize import fsolve
    x_min = -2.0 * np.pi
    x_max =  2.0 * np.pi
    y_min = -1.0
    y_max =  1.0
    z_min = -2.0 / 3.0 * np.pi
    z_max =  2.0 / 3.0 * np.pi
    dy_min = 0.016
    num_uniform = 7
    g = p3d.Grid().set_size(size, True)
    x = np.linspace(x_min, x_max, g.size[0,0] + 1)[:-1]
    z = np.linspace(z_min, z_max, g.size[0,2] + 1)[:-1]

    xi = np.linspace(-1.0, 1.0, g.size[0,1])
    y = np.sin(eta*xi*np.pi/2.0) / np.sin(eta*np.pi/2.0)
    y *= (y_max - y_min)/2.0

    dx, dz = (x_max-x_min)/g.size[0,0], (z_max-z_min)/g.size[0,2]
    dy = y[1:] - y[:-1]
    print ('dx : %.5E, dz : %.5E, dy_min : %.5E, dy_max : %.5E' % (dx,dz,np.min(dy),np.max(dy)) )

    for i in range(x.size):
        g.xyz[0][i,:,:,0] = x[i]
    for j in range(y.size):
        g.xyz[0][:,j,:,1] = y[j]
    for k in range(z.size):
        g.xyz[0][:,:,k,2] = z[k]
    return g

def initial_condition(g, mach_number=1.5, gamma=1.4):
    y = g.xyz[0][:,:,:,1]
    u0 = mach_number * 1.5 * ( 1.0 - y**2 )

    s = p3d.Solution().copy_from(g).quiescent(gamma)
    s.q[0][:,:,:,0] = 1.0
    s.q[0][:,:,:,2:4] = 0.0
    s.q[0][:,:,:,4] = 1.0 / gamma

    s.q[0][:,:,:,1] = u0
    for k in range(7):
        amp = 0.1 * np.random.random()
        s.q[0][:,1:-1,:,1] += u0[:,1:-1,:] * amp * np.sin( (2*k+3) * np.pi * y[:,1:-1,:] )
    return s.fromprimitive()

if __name__ == '__main__':
    g = grid([300, 200, 300],eta=0.99)
    g.save('ChannelFlow.xyz')
    initial_condition(g).save('ChannelFlow.ic.q')
