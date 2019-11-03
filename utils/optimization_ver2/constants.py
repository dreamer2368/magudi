NumSearch = 10
Nsplit = 3                                               # number of split time segments
Nts = 6                                                  # number of timesteps of each time segment
NtimestepOffset = 4                                      # timestep offset between time segments
startTimestep = 12                                       # initial timestep of the first time segment
totalTimestep = (Nsplit-1) * NtimestepOffset + Nts       # number of timesteps for the entire time span

matchingConditionWeight = 1.0e-7                         # weight for matching condition penalty
initialConditionControllability = 1.0e0                  # weight for derivative with respect to initial conditions

NcontrolRegion = 1                                       # number of control region

NcontrolSpace = NcontrolRegion + Nsplit - 1              # dimension of control space

initial_step = 1.0e9
golden_ratio = 1.618034
tol, eps = 1.0e-2,  1.0e-15
