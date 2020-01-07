from base import *

zeroControlForcing = True

command = ''
command += forwardRunCommand('.',zeroControlForcing)
command += adjointRunCommand()
command += gatherControlForcingGradientCommand()
command += innerProductCommand(globalGradFiles,globalGradFiles,ggFiles)
command += 'python3 optimization.py 1 -zero_baseline -initial_cg\n'
#command += 'python3 optimization.py 1 -initial_cg\n'

fID = open('initial.sh','w')
fID.write(command)
fID.close()

