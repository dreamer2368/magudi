Nsplit = 3                                               # number of split time segments
Nts = 6                                                  # number of timesteps of each time segment
startTimestep = 12                                       # initial timestep of the first time segment
totalTimestep = Nsplit * Nts                             # number of timesteps for the entire time span

import numpy as np
matchingConditionWeight = 1.0e-7 * np.ones(Nsplit)       # weight for matching condition penalty
initialConditionControllability = 1.0e0 * np.ones(Nsplit)# weight for derivative with respect to initial conditions
useLagrangian = False                                    # flag for augmented lagrangian
ignoreObjective = False
periodicSolution = False
if (not periodicSolution):
    matchingConditionWeight[0] = 0.0
    initialConditionControllability[0] = 0.0

NcontrolRegion = 1                                       # number of control region

NcontrolSpace = NcontrolRegion + Nsplit                  # dimension of control space

scriptorType = 'flux'                                    # either 'bash' or 'flux'
enableParallelBash = True
bashVerbose = False                                      # check only when serial bash loops
pcc = 36
maxNodes = 10

procedureSwitcher = { # number of nodes and processors for each procedure.
    'forward':      3 * np.array([1,pcc]),
    'adjoint':      3 * np.array([1,pcc]),
    'zaxpy':        3 * np.array([1,pcc]),
    'qfile-zaxpy':  3 * np.array([1,pcc]),
    'zxdoty':       3 * np.array([1,pcc]),
    'qfile-zxdoty': 3 * np.array([1,pcc]),
    'paste':        3 * np.array([1,pcc]),
    'slice':        3 * np.array([1,pcc]),
}

initial_step = 1.0e9
golden_ratio = 1.618034
tol, eps = 1.0e-2,  1.0e-15
linminTol, Nlinmin = 1.0e-1, 50
