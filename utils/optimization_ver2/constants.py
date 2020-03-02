Nsplit = 3                                               # number of split time segments
Nts = 6                                                  # number of timesteps of each time segment
NtimestepOffset = 4                                      # timestep offset between time segments
startTimestep = 12                                       # initial timestep of the first time segment
totalTimestep = (Nsplit-1) * NtimestepOffset + Nts       # number of timesteps for the entire time span

import numpy as np
matchingConditionWeight = 1.0e-7 * np.ones(Nsplit)       # weight for matching condition penalty
matchingConditionWeight[-1] = 0.0
initialConditionControllability = 1.0e0 * np.ones(Nsplit)# weight for derivative with respect to initial conditions
initialConditionControllability[0] = 0.0

NcontrolRegion = 1                                       # number of control region

NcontrolSpace = NcontrolRegion + Nsplit                  # dimension of control space

enableParallelBash = True
bashVerbose = False                                      # check only when serial bash loops
pcc = 36
maxNodes = 10

NprocForward, NprocAdjoint = pcc*3, pcc*3                # number of processors for forward/adjoint run
NprocZaxpy, NprocQfileZaxpy = 10, 10                     # number of processors for zaxpy works
NprocZxdoty, NprocQfileZxdoty = 10, 10                   # number of processors for inner product works
NprocPaste, NprocSlice = 10, 10                          # number of processors for paste/slice control forcing

NodesForward, NodesAdjoint = 3, 3                        # number of nodes for forward/adjoint run
NodesZaxpy, NodesQfileZaxpy = 1, 1                       # number of nodes for zaxpy works
NodesZxdoty, NodesQfileZxdoty = 1, 1                     # number of nodes for inner product works
NodesPaste, NodesSlice = 1, 1                            # number of nodes for paste/slice control forcing

initial_step = 1.0e9
golden_ratio = 1.618034
tol, eps = 1.0e-2,  1.0e-15
linminTol = 1.0e-1
