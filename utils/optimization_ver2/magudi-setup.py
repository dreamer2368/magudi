from base import *

commandFile = 'magudi-setup.sh'
commandString = ''
commandString += generalSetOptionCommand

commands = []
for bdir in ['a','b','c','x','x0']:
    targetInputFiles = ['%s/%s/%s'%(bdir,dir,file) for dir, file in zip(directories,inputFiles)]
    for k in range(Nsplit):
        commands += ['./setOption.sh %s "output_prefix" "\'WeiFreundSDML-%d\'" ' % (targetInputFiles[k],k) ]
commandString += scriptor.nonMPILoopCommand(commands,'output_prefix')

commands = []
for bdir in ['a','b','c','x','x0']:
    targetInputFiles = ['%s/%s/%s'%(bdir,dir,file) for dir, file in zip(directories,inputFiles)]
    for k in range(Nsplit):
        commands += ['./setOption.sh %s "initial_condition_file" "\'WeiFreundSDML-%d.ic.q\'" ' % (targetInputFiles[k],k)]
commandString += scriptor.nonMPILoopCommand(commands,'initial_condition')

commands = []
for bdir in ['x0']:
    targetInputFiles = ['%s/%s/%s'%(bdir,dir,file) for dir, file in zip(directories,inputFiles)]
    for k in range(Nsplit):
        commands += ['./setOption.sh %s "controller_switch" "false" ' % targetInputFiles[k]]
commandString += scriptor.nonMPILoopCommand(commands,'controller_switch_false')

#commands = []
#for bdir in ['a','b','c','x','x0']:
#    targetInputFiles = ['%s/%s/%s'%(bdir,dir,file) for dir, file in zip(directories,inputFiles)]
#    for k in range(Nsplit):
#        commandString += './setOption.sh %s "output_prefix" "\'WeiFreundSDML-%d\'" \n'         \
#                        % (targetInputFiles[k],k)
#        commandString += './setOption.sh %s "initial_condition_file" "\'..\/WeiFreundSDML-%d.ic.q\'" \n'         \
#                        % (targetInputFiles[k],k)
#        #commandString += './setOption.sh %s "controller_switch" "true" \n' % targetInputFiles[k]
#        #commandString += './setOption.sh %s "controller_factor" "1.0e-3" \n' % (targetInputFiles[k])
#        #commandString += './setOption.sh %s "controller_norm" "\'L_Inf_with_timestep\'" \n' % (targetInputFiles[k])

#for bdir in ['x0']:
#    targetInputFiles = ['%s/%s/%s'%(bdir,dir,file) for dir, file in zip(directories,inputFiles)]
#    for k in range(Nsplit):
#        commandString += './setOption.sh %s "controller_switch" "false" \n' % targetInputFiles[k]
#        #commandString += './setOption.sh %s "report_interval" 50 \n' % targetInputFiles[k]
#        #commandString += './setOption.sh %s "save_interval" 100 \n' % targetInputFiles[k]
#        #commandString += './setOption.sh %s "controller_buffer_size" 400 \n' % targetInputFiles[k]
#        #commandString += './setOption.sh %s "adjoint_nonzero_initial_condition" false \n' % targetInputFiles[k]

fID = open(commandFile,'w')
fID.write(commandString)
fID.close()
