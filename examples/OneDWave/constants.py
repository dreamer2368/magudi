Nsplit = 24                                              # number of split time segments
Nts = 20                                                 # number of timesteps of each time segment
startTimestep = 0                                        # initial timestep of the first time segment
totalTimestep = Nsplit * Nts                             # number of timesteps for the entire time span

import numpy as np
penaltyType = 'Huber'                                    # either 'base' or 'Huber'
matchingConditionWeight = 1.0e-5 * np.ones(Nsplit)       # weight for matching condition penalty
initialConditionControllability = 1.0e0 * np.ones(Nsplit)# weight for derivative with respect to initial conditions
useLagrangian = False                                    # flag for augmented lagrangian
ignoreIntegralObjective = False
terminalObjective = False                                # flag for terminal objective
ignoreController = False
periodicSolution = False
if (not periodicSolution):
    matchingConditionWeight[0] = 0.0
    initialConditionControllability[0] = 0.0

NcontrolRegion = 1                                       # number of control region
if (ignoreController):
    NcontrolRegion = 0

NcontrolSpace = NcontrolRegion + Nsplit                  # dimension of control space

scriptorType = 'base'                                    # 'base', 'bash' or 'flux'
enableParallelBash = True
bashVerbose = True if (enableParallelBash) else False    # check only when serial bash loops
pcc = 1
maxNodes = 1
maxProcs = maxNodes * pcc

procedureSwitcher = { # number of nodes and processors for each procedure.
    'forward':      np.array([1,pcc]),
    'adjoint':      np.array([1,pcc]),
    'zaxpy':        np.array([1,pcc]),
    'qfile-zaxpy':  np.array([1,1]),
    'zxdoty':       np.array([1,pcc]),
    'qfile-zxdoty': np.array([1,1]),
    'paste':        np.array([1,1]),
    'slice':        np.array([1,1]),
}

initial_step, safe_zone = 1.0e0, 1.0e4
golden_ratio = 1.618034
tol, eps, huge = 1.0e-9,  1.0e-15, 1.0e100
linminTol, Nlinmin = 1.0e-1, 50

saveDiffFiles = False                                     # flag for saving discontinuity files
Ndiff = 8                                                # number of discontinuity files that are saved in line searches
diffStep = [1,2,3,4,5,6,10,15]

saveStateLog = False
NdiffRegion = 4
