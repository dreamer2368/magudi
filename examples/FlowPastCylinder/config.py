#!/usr/bin/env python
import numpy as np
import plot3dnasa as p3d

def grid(size):
    x_min =   0.
    x_max =  30.
    y_min =  -5.
    y_max =   5.
    g = p3d.Grid().set_size(size, True)
    x = np.linspace(x_min, x_max, g.size[0,0])
    y = np.linspace(y_min, y_max, g.size[0,1])
    g.xyz[0][:,:,0,:2] = np.transpose(np.meshgrid(x, y))
    return g

def initial_condition(g, Ma, gamma=1.4):
    s = p3d.Solution().copy_from(g).quiescent(gamma)

    s.q[0][:,:,:,1] = Ma
    return s.fromprimitive(gamma)

if __name__ == '__main__':
    g = grid([3011, 1011])
    g.save('FlowPastCylinder.xyz')
    Ma, gamma = 0.1, 1.4
    s = initial_condition(g, Ma, gamma)
    s.save('FlowPastCylinder.ic.q')
    s.save('FlowPastCylinder.target.q')

