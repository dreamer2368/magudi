from constants import *
from filenames import *

import numpy as np
import subprocess
import pandas as pd

def bashCheckResultCommand(procedureName, countFails=0):
    if(countFails>0):
        commandString = 'FAIL=0\n'
        commandString += 'for k in {0..%d}\n'%(countFails-1)
        commandString += 'do\n'                                                             \
                         '    wait ${pids[${k}]} || let "FAIL+=1"\n'                        \
                         '    echo ${pids[${k}]}\n'                                         \
                         'done\n\n'
        commandString += 'echo "Number of failures: $FAIL"\n'                                \
                         'if [ $FAIL -ne 0 ]; then\n'
    else:
        commandString = 'if [ $? -ne 0 ]; then\n'

    commandString += '   echo "'+procedureName+' failed."\n'                                \
                     '   exit -1\n'                                                         \
                     'else\n'                                                               \
                     '   echo "'+procedureName+' succeeded."\n'
    if(countFails>0):
        commandString += '   FAIL=0\n'
    commandString += 'fi\n\n'

    return commandString

def bashGetNodeListCommand():
    commandString = 'nodeList=()\n'                                                         \
                    'while IFS= read -r line\n'                                             \
                    'do\n'                                                                  \
                    '    nodeList+=("$line")\n'                                             \
                    'done < <( scontrol show hostnames ${SLURM_JOB_NODELIST} )\n'
    return commandString

def bashGetNodeListSliceCommand(index, numNodes):
    nodeIndex = numNodes * index
    commandString = 'let "nodeIndex=%d"\n' % nodeIndex
    commandString += 'nodeListString="${nodeList[${nodeIndex}]}"\n'
    if (numNodes>1):
        commandString += 'for j in {1..%d}\n' % (numNodes-1)
        commandString += 'do\n'
        commandString += '    let "nodeIndex=%d+${j}"\n' % nodeIndex
        commandString += '    nodeListString+=",${nodeList[${nodeIndex}]}"\n'               \
                         'done\n\n'

    return commandString
