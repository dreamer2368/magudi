from base import *

zeroControlForcing = True
initialLagrangian = True
weight = 0.0e-1 * np.ones(Nsplit)       # weight for matching condition penalty
if (not periodicSolution):
    weight[0] = 0.0

command = ''
command += forwardRunCommand('x0',zeroControlForcing)
command += updateLagrangian(weight,initialLagrangian)

fID = open('initial-forward.sh','w')
fID.write(command)
fID.close()

command = ''
command += adjointRunCommand()
command += gatherControlForcingGradientCommand()
command += innerProductCommand(globalGradFiles,globalGradFiles,ggFiles)

fID = open('initial-adjoint.sh','w')
fID.write(command)
fID.close()
