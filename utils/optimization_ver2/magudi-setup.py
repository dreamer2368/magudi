from base import *

commandFile = 'magudi-setup.sh'
commandString = ''
commandString += generalSetOptionCommand

for bdir in ['x0']:
    targetInputFiles = ['%s/%s/%s'%(bdir,dir,file) for dir, file in zip(directories,inputFiles)]
    for k in range(21):
        commandString += 'setOption %s "initial_condition_file" "WeiFreundSDML-%d.ic.q" \n'         \
                        % (targetInputFiles[k],k)

for bdir in ['a','b','c','x']:
    targetInputFiles = ['%s/%s/%s'%(bdir,dir,file) for dir, file in zip(directories,inputFiles)]
    for k in range(21):
        commandString += 'setOption %s "controller_switch" "true" \n' % targetInputFiles[k]

fID = open(commandFile,'w')
fID.write(commandString)
fID.close()
