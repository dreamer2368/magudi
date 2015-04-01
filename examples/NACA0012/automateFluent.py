#!/usr/bin/env python
import numpy as np

if __name__ == '__main__':

    outputPrefix = 'NACA0012'
    caseFile = outputPrefix + '.cas'

    nIterations = 10000

    freeStreamDensity = 1.225
    specificHeatAtConstantPressure = 1006.43
    gasConstant = 8.314462175 * 1000. / 28.966
    specificHeatAtConstantVolume = specificHeatAtConstantPressure - gasConstant
    ratioOfSpecificHeats = specificHeatAtConstantPressure / specificHeatAtConstantVolume
    freeStreamPressure = 101325.
    freeStreamSpeedOfSound = np.sqrt(ratioOfSpecificHeats * freeStreamPressure / freeStreamDensity)

    journalFile = 'fluentScript.jou'
    f = open(journalFile, 'w')

    print >>f, 'rc "%s"' % (caseFile)

    print >>f, """
file/set-batch-options no no yes no
grid/modify-zones/zone-name 3 interior
grid/modify-zones/zone-name 4 freestream
grid/modify-zones/zone-name 5 airfoil
define/models/viscous/inviscid yes
define/models/energy yes
"""

    print >>f, 'define/boundary-conditions/velocity-inlet freestream yes yes no %f no 1 no 0 no %f' % (0.2 * freeStreamSpeedOfSound, freeStreamPressure / (freeStreamDensity * gasConstant))

    print >>f, """
solve/initialize/compute-defaults/velocity-inlet freestream
solve/monitors/residual/convergence-criteria 1e-4 1e-7 1e-7 1e-9
solve/set/discretization-scheme/pressure 12
solve/set/discretization-scheme/mom 1
solve/set/discretization-scheme/enthalpy 1
solve/initialize/initialize-flow
solve/initialize/initialize-flow
"""

    print >>f, 'solve/iterate %i' % (nIterations)
    print >>f, 'wcd "%s"' % (outputPrefix)

    print >>f, 'file/export/ascii "%s" interior airfoil freestream () no pressure y-velocity x-velocity density q yes' % (outputPrefix + '_steady_ascii.txt')
    print >>f, 'exit'

    f.close()
