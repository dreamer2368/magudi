from base import *

commandFile = 'magudi-setup.sh'
commandString = ''

for k in range(21):
    magudiFile = '%d/magudi-%d.inp'%(k,k)
    commandString += setOptionCommand(magudiFile)

    commandString += 'setOption "initial_condition_file" "WeiFreundSDML-%d.ic.q" \n' % k

for dir in ['a','b','c','x']:
    for k in range(21):
        magudiFile = '%s/%d/magudi-%d.inp'%(dir,k,k)
        commandString += setOptionCommand(magudiFile)

        commandString += 'setOption "controller_switch" "true" \n'


fID = open(commandFile,'w')
fID.write(commandString)
fID.close()

