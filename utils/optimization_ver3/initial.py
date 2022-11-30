from inputs import InputParser
from optimizer import Optimizer
import numpy as np

config = InputParser('./optim.yml')
optim = Optimizer(config)

zeroControlForcing = True
initialLagrangian = True
weight = 0.0e-1 * np.ones(optim.const.Nsplit)       # weight for matching condition penalty
if (not optim.const.periodicSolution):
    weight[0] = 0.0

command = ''
command += optim.base.forwardRunCommand('x0',zeroControlForcing)
command += optim.base.updateLagrangian(weight,initialLagrangian)

fID = open('initial-forward.sh','w')
fID.write(command)
fID.close()

command = ''
command += optim.base.adjointRunCommand()
command += optim.base.gatherControlForcingGradientCommand()
command += optim.base.innerProductCommand(optim.fl.globalGradFiles,optim.fl.globalGradFiles,optim.fl.ggFiles)

fID = open('initial-adjoint.sh','w')
fID.write(command)
fID.close()
