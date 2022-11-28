from base import *

commandFile = 'magudi-setup.sh'
commandString = ''
commandString += generalSetOptionCommand

commands = []
for bdir in ['x','x0']:
    targetInputFiles = ['%s/%s/%s'%(bdir,dir,file) for dir, file in zip(directories,inputFiles)]
    for k in range(Nsplit):
        commands += ['./setOption.sh %s "output_prefix" "\'OneDWave-%d\'" ' % (targetInputFiles[k],k) ]
commandString += scriptor.nonMPILoopCommand(commands,'output_prefix')

commands = []
for bdir in ['x','x0']:
    targetInputFiles = ['%s/%s/%s'%(bdir,dir,file) for dir, file in zip(directories,inputFiles)]
    for k in range(Nsplit):
        commands += ['./setOption.sh %s "number_of_timesteps" "20"' % targetInputFiles[k]]
commandString += scriptor.nonMPILoopCommand(commands,'initial_condition')

commands = []
for bdir in ['x','x0']:
    targetInputFiles = ['%s/%s/%s'%(bdir,dir,file) for dir, file in zip(directories,inputFiles)]
    for k in range(Nsplit):
        commands += ['./setOption.sh %s "initial_condition_file" "\'OneDWave-%d.ic.q\'" ' % (targetInputFiles[k],k)]
commandString += scriptor.nonMPILoopCommand(commands,'initial_condition')

commands = []
for bdir in ['x']:
    targetInputFiles = ['%s/%s/%s'%(bdir,dir,file) for dir, file in zip(directories,inputFiles)]
    for k in range(Nsplit):
        commands += ['./setOption.sh %s "controller_switch" "true" ' % targetInputFiles[k]]
commandString += scriptor.nonMPILoopCommand(commands,'controller_switch_true')

fID = open(commandFile,'w')
fID.write(commandString)
fID.close()
