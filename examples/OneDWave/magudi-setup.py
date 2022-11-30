from inputs import InputParser
from optimizer import Optimizer

config = InputParser('./optim.yml')
optim = Optimizer(config)

commandFile = 'magudi-setup.sh'
commandString = ''
commandString += optim.base.generalSetOptionCommand

commands = []
for bdir in ['a', 'b', 'c', 'x','x0']:
    targetInputFiles = ['%s/%s/%s'%(bdir,dir,file) for dir, file in zip(optim.fl.directories, optim.fl.inputFiles)]
    for k in range(optim.const.Nsplit):
        commands += ['./setOption.sh %s "output_prefix" "\'OneDWave-%d\'" ' % (targetInputFiles[k],k) ]
commandString += optim.scriptor.nonMPILoopCommand(commands,'output_prefix')

commands = []
for bdir in ['a', 'b', 'c', 'x','x0']:
    targetInputFiles = ['%s/%s/%s'%(bdir,dir,file) for dir, file in zip(optim.fl.directories, optim.fl.inputFiles)]
    for k in range(optim.const.Nsplit):
        commands += ['./setOption.sh %s "number_of_timesteps" "20"' % targetInputFiles[k]]
commandString += optim.scriptor.nonMPILoopCommand(commands,'initial_condition')

commands = []
for bdir in ['a', 'b', 'c', 'x','x0']:
    targetInputFiles = ['%s/%s/%s'%(bdir,dir,file) for dir, file in zip(optim.fl.directories, optim.fl.inputFiles)]
    for k in range(optim.const.Nsplit):
        commands += ['./setOption.sh %s "initial_condition_file" "\'OneDWave-%d.ic.q\'" ' % (targetInputFiles[k],k)]
commandString += optim.scriptor.nonMPILoopCommand(commands,'initial_condition')

commands = []
for bdir in ['a', 'b', 'c', 'x']:
    targetInputFiles = ['%s/%s/%s'%(bdir,dir,file) for dir, file in zip(optim.fl.directories, optim.fl.inputFiles)]
    for k in range(optim.const.Nsplit):
        commands += ['./setOption.sh %s "controller_switch" "true" ' % targetInputFiles[k]]
commandString += optim.scriptor.nonMPILoopCommand(commands,'controller_switch_true')

fID = open(commandFile,'w')
fID.write(commandString)
fID.close()
