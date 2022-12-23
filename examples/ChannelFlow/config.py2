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
    x = g.xyz[0][:,:,:,0]
    y = g.xyz[0][:,:,:,1]
    z = g.xyz[0][:,:,:,2]
    u0 = mach_number * 1.5 * ( 1.0 - y**2 )

    s = p3d.Solution().copy_from(g).quiescent(gamma)
    s.q[0][:,:,:,0] = 1.0
    s.q[0][:,:,:,1:4] = 0.0
    s.q[0][:,:,:,4] = 1.0 / gamma

    temp = np.zeros_like(s.q[0][:,:,:,1:4],dtype=np.complex128)
    for k in range(3,6):
        for l in range(18,21):
            for m in range(6,9):
                for dim in range(3):
                    amp = np.random.rand(2)
                    amp = amp[0] + amp[1] * 1.0j
                    phase = 2.0 * np.pi * np.random.rand(2)
                    yphase = 2 * np.random.randint(2) - 1
                    print ('amp: (%.15E + %.15E i), phase: (%.15E, %.15E), y-phase: %d'
                            %(amp.real,amp.imag,phase[0],phase[1],yphase))
                    temp[:,1:-1,:,dim] += u0[:,1:-1,:] * amp                    \
                              * np.sin( k * np.pi * y[:,1:-1,:] * yphase )      \
                      * np.exp( l * 0.5j * x[:,1:-1,:] + phase[0] * 1.0j )      \
                      * np.exp( m * 1.5j * z[:,1:-1,:] + phase[1] * 1.0j )

    temp = temp.real
    temp *= 0.05 * mach_number * 1.5 / np.amax(np.abs(temp))

    s.q[0][:,:,:,1] = u0
    s.q[0][:,:,:,1:4] += temp.real
    return s.fromprimitive()

if __name__ == '__main__':
    g = grid([300, 200, 300],eta=0.99)
    g.save('ChannelFlow.xyz')
    initial_condition(g).save('ChannelFlow.ic.q')
