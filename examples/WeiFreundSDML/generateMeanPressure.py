#!/usr/bin/env python
import numpy as np

import PLOT3D

if __name__ == '__main__':

    outputPrefix = 'WeiFreundSDML'

    ratioOfSpecificHeats = 1.4

    soln = PLOT3D.Solution()
    soln.Import(outputPrefix + '.mean.q')

    func = PLOT3D.Function()
    func.nComponents = 1
    func.CopyFrom(soln)

    func.F[0][:,:,:,0] = (ratioOfSpecificHeats - 1.) * (soln.Q[0][:,:,:,-1] - 0.5 / soln.Q[0][:,:,:,0] * np.sum(soln.Q[0][:,:,:,1:-1] ** 2, axis = -1))
    func.Export(outputPrefix + '.mean_pressure.f')
