#!/usr/bin/env python
import numpy as np
from scipy.optimize import fsolve
import numpy.random as random

import plot3dnasa as p3d

EPSILON = np.finfo(float).eps

def grid():
    Lx, Ly, Lz = 10 * random.rand(3,1)
    x_min = -Lx
    x_max =  Lx
    y_min = -Ly
    y_max =  Ly
    z_min = -Lz
    z_max =  Lz

    numGrid = random.randint(20, size=[3,]) + 20
    numX, numY, numZ = numGrid[0], numGrid[1], numGrid[2]
    size = [[numX,numY,numZ],
            [numY,numX,numZ],
            [numX,numY,numZ],
            [numY,numX,numZ],
            [numX,numY,numZ]]
    g = p3d.Grid().set_size(size,True)
    print (g.size)

    dx = [ 2.0*Lx/(numX-1), 2.0*Ly/(numY-1), 2.0*Lz/(numZ-1) ]

    x0 = np.linspace(x_min,x_max, g.size[0,0])
    y0 = np.linspace(y_min,y_max, g.size[0,1])
    z0 = np.linspace(z_min,z_max, g.size[0,2])

    # center, N, W, S, E
    g.xyz[0] = np.moveaxis(np.meshgrid(x0, y0, z0, indexing='ij'       ),[0,1,2,3],[3,0,1,2])
    g.xyz[1] = np.moveaxis(np.meshgrid(-x0, y0+2*Ly, z0, indexing='ij' ),[0,1,2,3],[3,1,0,2])
    g.xyz[2] = np.moveaxis(np.meshgrid(-2*Lx-x0, -y0, z0, indexing='ij'),[0,1,2,3],[3,0,1,2])
    g.xyz[3] = np.moveaxis(np.meshgrid(x0, -y0-2*Ly, z0, indexing='ij' ),[0,1,2,3],[3,1,0,2])
    g.xyz[4] = np.moveaxis(np.meshgrid(2*Lx+x0, y0, z0, indexing='ij'  ),[0,1,2,3],[3,0,1,2])

    gridPerturbation = random.rand(*np.shape(g.xyz[0]))
    gridPerturbation = 2.0 * gridPerturbation - 1.0
    for k in range(3):
        gridPerturbation[:,:,:,k] *= 0.1 * dx[k]

    g.xyz[0] += gridPerturbation
    g.xyz[1][0,:,:,:] += gridPerturbation[::-1,-1,:,:]
    g.xyz[2][0,:,:,:] += gridPerturbation[0,::-1,:,:]
    g.xyz[3][0,:,:,:] += gridPerturbation[:,0,:,:]
    g.xyz[4][0,:,:,:] += gridPerturbation[-1,:,:,:]

    for i in range(1,5):
        gridPerturbation = random.rand(*np.shape(g.xyz[i]))
        gridPerturbation = 2.0 * gridPerturbation - 1.0
        for k in range(3):
            gridPerturbation[:,:,:,k] *= 0.1 * dx[k]
        gridPerturbation[0,:,:,:] *= 0.0
        g.xyz[i] += gridPerturbation

    idxX = random.randint(numX)
    idxY = random.randint(numY)
    idxZ = random.randint(numZ)
    print ('0-1 interface')
    print (g.xyz[0][idxX,-1,idxZ,:])
    print (g.xyz[1][0,numX-idxX-1,idxZ,:])

    print ('0-2 interface')
    print (g.xyz[0][0,idxY,idxZ,:])
    print (g.xyz[2][0,numY-idxY-1,idxZ,:])

    print ('0-3 interface')
    print (g.xyz[0][idxX,0,idxZ,:])
    print (g.xyz[3][0,idxX,idxZ,:])

    print ('0-4 interface')
    print (g.xyz[0][-1,idxY,idxZ,:])
    print (g.xyz[4][0,idxY,idxZ,:])

    print ('center grid test points indexes')
    print (np.array([[idxX,numY-1,idxZ],
                     [0,idxY,idxZ],
                     [idxX,0,idxZ],
                     [numX-1,idxY,idxZ]])+1)

    print ('patches/interface.N1/test_index_1 = %d'%(idxX+1))
    print ('patches/interface.N1/test_index_2 = %d'%(numY))
    print ('patches/interface.N1/test_index_3 = %d'%(idxZ+1))

    print ('patches/interface.W1/test_index_1 = %d'%(1))
    print ('patches/interface.W1/test_index_2 = %d'%(idxY+1))
    print ('patches/interface.W1/test_index_3 = %d'%(idxZ+1))

    print ('patches/interface.S1/test_index_1 = %d'%(idxX+1))
    print ('patches/interface.S1/test_index_2 = %d'%(1))
    print ('patches/interface.S1/test_index_3 = %d'%(idxZ+1))

    print ('patches/interface.E1/test_index_1 = %d'%(numX))
    print ('patches/interface.E1/test_index_2 = %d'%(idxY+1))
    print ('patches/interface.E1/test_index_3 = %d'%(idxZ+1))


    print ('circumferential grid test points indexes')
    print (np.array([[0,numX-idxX-1,idxZ],
                     [0,numY-idxY-1,idxZ],
                     [0,idxX,idxZ],
                     [0,idxY,idxZ]])+1)

    print ('patches/interface.N2/test_index_1 = %d'%(1))
    print ('patches/interface.N2/test_index_2 = %d'%(numX-idxX))
    print ('patches/interface.N2/test_index_3 = %d'%(idxZ+1))

    print ('patches/interface.W3/test_index_1 = %d'%(1))
    print ('patches/interface.W3/test_index_2 = %d'%(numY-idxY))
    print ('patches/interface.W3/test_index_3 = %d'%(idxZ+1))

    print ('patches/interface.S4/test_index_1 = %d'%(1))
    print ('patches/interface.S4/test_index_2 = %d'%(idxX+1))
    print ('patches/interface.S4/test_index_3 = %d'%(idxZ+1))

    print ('patches/interface.E5/test_index_1 = %d'%(1))
    print ('patches/interface.E5/test_index_2 = %d'%(idxY+1))
    print ('patches/interface.E5/test_index_3 = %d'%(idxZ+1))

    return g


if __name__ == '__main__':
    g = grid()
    output_prefix = 'FiveBlockMesh'
    g.save(output_prefix + '.xyz')
