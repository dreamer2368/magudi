class Constants:
    Nsplit = -1                 # number of split time segments
    Nts = -1                    # number of timesteps of each time segment
    startTimestep = -1          # initial timestep of the first time segment
    totalTimestep = -1          # number of timesteps for the entire time span

    def __init__(self, config):
        self.Nsplit = config.getInput(['time_splitting', 'number_of_segments'], datatype=int)
        self.Nts = config.getInput(['time_splitting', 'segment_length'], datatype=int)
        self.startTimestep = config.getInput(['time_splitting', 'start_timestep'], datatype=int)
        self.totalTimestep = self.Nsplit * self.Nts
        return
# Nsplit = 6                                               # number of split time segments
# Nts = 2400                                               # number of timesteps of each time segment
# startTimestep = 30000                                    # initial timestep of the first time segment
# totalTimestep = Nsplit * Nts                             # number of timesteps for the entire time span

import numpy as np
# penaltyType = 'base'                                     # either 'base' or 'Huber'
matchingConditionWeight = 1.6e-6 * np.ones(Nsplit)       # weight for matching condition penalty
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

# scriptorType = 'flux'                                    # 'base', 'bash' or 'flux'
# enableParallelBash = True
# bashVerbose = False                                      # check only when serial bash loops
# pcc = 36
# maxNodes = 6
# maxProcs = maxNodes * pcc

# procedureSwitcher = { # number of nodes and processors for each procedure.
#     'forward':      maxNodes * np.array([1,pcc]),
#     'adjoint':      maxNodes * np.array([1,pcc]),
#     'zaxpy':        2 * np.array([1,pcc]),
#     'qfile-zaxpy':  1 * np.array([1,5]),
#     'zxdoty':       2 * np.array([1,pcc]),
#     'qfile-zxdoty': 1 * np.array([1,5]),
#     'paste':        1 * np.array([1,5]),
#     'slice':        1 * np.array([1,5]),
# }

initial_step, safe_zone = 1.0e0, 1.0e4
golden_ratio = 1.618034
tol, eps, huge = 1.0e-8,  1.0e-15, 1.0e100
linminTol, Nlinmin = 1.0e-1, 50

saveDiffFiles = False                                    # flag for saving discontinuity files
Ndiff = 8                                                # number of discontinuity files that are saved in line searches
diffStep = int(np.floor(Nsplit/Ndiff))

saveStateLog = True
NdiffRegion = 4
