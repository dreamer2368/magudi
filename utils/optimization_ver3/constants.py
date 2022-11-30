import numpy as np

class Constants:
    Nsplit = -1                 # number of split time segments
    Nts = -1                    # number of timesteps of each time segment
    startTimestep = -1          # initial timestep of the first time segment
    totalTimestep = -1          # number of timesteps for the entire time span

    ignoreController = False
    NcontrolRegion = -1         # number of control region
    NcontrolSpace = -1          # dimension of control space including intermediate states

    saveStateLog = True
    NdiffRegion = -1

    ignoreIntegralObjective = False

    penaltyType = ''
    useLagrangian = False       # flag for augmented lagrangian
    terminalObjective = False   # flag for terminal objective

    periodicSolution = False
    matchingConditionWeight = 1.6e-6 * np.ones(Nsplit)       # weight for matching condition penalty
    initialConditionControllability = 1.0e0 * np.ones(Nsplit)# weight for derivative with respect to initial conditions

    saveDiffFiles = False                                    # flag for saving discontinuity files
    Ndiff = -1                                                # number of discontinuity files that are saved in line searches
    diffStep = -1

    initial_step, safe_zone = -1.0, -1.0
    tol = -1.0
    linminTol, Nlinmin = -1.0, -1

    golden_ratio, eps, huge = 1.618034, 1.0e-15, 1.0e100

    def __init__(self, config):
        self.Nsplit = config.getInput(['time_splitting', 'number_of_segments'], datatype=int)
        self.Nts = config.getInput(['time_splitting', 'segment_length'], datatype=int)
        self.startTimestep = config.getInput(['time_splitting', 'start_timestep'], datatype=int)
        self.totalTimestep = self.Nsplit * self.Nts

        self.ignoreController = (not config.getInput(['controller', 'enabled'], datatype=bool))
        if (self.ignoreController):
            self.NcontrolRegion = 0
        else:
            self.NcontrolRegion = config.getInput(['controller', 'number_of_actuators', datatype=int])
            assert(self.NcontrolRegion > 0)

        self.NcontrolSpace = self.NcontrolRegion + self.Nsplit

        self.saveStateLog = config.getInput(['state_log', 'enabled'], fallback=False)
        self.NdiffRegion = config.getInput(['state_log', 'number_of_regions'], fallback=0)

        self.ignoreIntegralObjective = (not config.getInput(['objective', 'include_time_integral'], fallback=True))
        self.terminalObjective = config.getInput(['objective', 'include_terminal'], fallback=False)

        self.penaltyType = config.getInput(['penalty_norm', 'base'], datatype=str)
        if (self.penaltyType = 'Huber'):
            self.useLagrangian = False
        else:
            self.useLagrangian = config.getInput(['penalty_norm', 'augmented_lagrangian'], fallback=False)

        self.matchingConditionWeight = np.ones(self.Nsplit)
        self.matchingConditionWeight *= config.getInput(['time_splitting', 'matching_condition_weight'], datatype=float)
        self.initialConditionControllability = np.ones(self.Nsplit)
        self.initialConditionControllability *= config.getInput(['time_splitting', 'state_controllability'], fallback=1.0)
        self.periodicSolution = config.getInput(['time_splitting', 'periodic_solution'], fallback=False)
        if (not self.periodicSolution):
            self.matchingConditionWeight[0] = 0.0
            self.initialConditionControllability[0] = 0.0

        self.saveDiffFiles = config.getInput(['discontinuity_log', 'enabled'], fallback=False)
        if (self.saveDiffFiles):
            self.Ndiff = config.getInput(['discontinuity_log', 'number_of_time_points'], datatype=int)
        else:
            self.Ndiff = 0
        self.diffStep = int(np.floor(self.Nsplit/self.Ndiff))

        self.initial_step = config.getInput(['optimization', 'line_minimization', 'initial_step_size'], fallback = 1.0)
        self.safe_zone = config.getInput(['optimization', 'line_minimization', 'safe_zone'], fallback = 1.0e4)
        self.tol = config.getInput(['optimization', 'tolerance'], fallback=1.0e-8)
        self.linminTol = config.getInput(['optimization', 'line_minimization', 'tolerance'], fallback=1.0e-1)
        self.Nlinmin = config.getInput(['optimization', 'line_minimization', 'number_of_searches'], fallback=50)
        return
