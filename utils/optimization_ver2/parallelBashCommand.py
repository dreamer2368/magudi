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
                         '    echo "${pids[${k}]}, FAIL: ${FAIL}"\n'                        \
                         'done\n\n'
        commandString += 'echo "Number of failures: $FAIL"\n'                               \
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
    commandString =  'if [ ${SLURM_JOB_NUM_NODES} -ne %d ]; then\n' % maxNodes
    commandString += '   echo "${SLURM_JOB_NUM_NODES} number of nodes are assigned."\n'     \
                     '   echo "max nodes in python: %d."\n' % maxNodes
    commandString += '   exit -1\n fi \n\n'
    commandString += 'nodeList=()\n'                                                        \
                     'while IFS= read -r line\n'                                            \
                     'do\n'                                                                 \
                     '    nodeList+=("$line")\n'                                            \
                     'done < <( scontrol show hostnames ${SLURM_JOB_NODELIST} )\n\n'
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

def bashSerialLoopCommand(commands,procedure,prefix='job',directories=None):
    if (directories is None):
        moveDir = False
    else:
        moveDir = True
        if (len(commands)!=len(directories)):
            raise ValueError('Provide directories for all commands.')

    nodePerCommand, procsPerCommand = procedureSwitcher.get(procedure)

    commandString = ''
    commandString += bashGetNodeListCommand()

    for k, command in enumerate(commands):
        if (moveDir): commandString += 'cd %s \n'%(directories[k])
        commandString += bashGetNodeListSliceCommand(0,nodePerCommand)
        commandString += 'srun -N %d -n %d -w ${nodeListString} '                         \
                         % (nodePerCommand,procsPerCommand)
        commandString += command
        commandString += ' &> %s/%s_result_%d.out \n' % (OUTDIR,prefix,k)
        commandString += bashCheckResultCommand('%s-serial-%d'%(prefix,k))
        if (moveDir): commandString += 'cd %s \n'%(ROOTDIR)

    return commandString

def bashParallelLoopCommand(commands,procedure,prefix='job',directories=None):
    if (directories is None):
        moveDir = False
    else:
        moveDir = True
        if (len(commands)!=len(directories)):
            raise ValueError('Provide directories for all commands.')

    nodePerCommand, procsPerCommand = procedureSwitcher.get(procedure)

    maxJobsPerLoop = int(np.floor(maxNodes/nodePerCommand))
    if (maxJobsPerLoop<1):
        raise ValueError('The job %s requires at least %d number of nodes. '                \
                            %(prefix,nodePerCommand) +                                      \
                         'Only %d nodes are assigned.'%(maxNodes))

    commandString = ''

    if (enableParallelBash): commandString += bashGetNodeListCommand()
    idx, nJobs = 0, len(commands)
    loop = 0
    while (nJobs>0):
        jobsPerLoop = maxJobsPerLoop if (nJobs>maxJobsPerLoop) else nJobs
        for k in range(jobsPerLoop):
            if (moveDir): commandString += 'cd %s \n'%(directories[idx+k])

            command = ''
            if (enableParallelBash):
                command += bashGetNodeListSliceCommand(k,nodePerCommand)
                command += 'srun -N %d -n %d -w ${nodeListString} '                         \
                            % (nodePerCommand,procsPerCommand)
            command += commands[idx+k]
            if (enableParallelBash or (not bashVerbose) ):
                command += ' &> %s/%s_result_%d.out &' % (OUTDIR,prefix,idx+k)
            command += '\n'
            if (enableParallelBash):
                command += 'pids[%d]=$!\n\n' % k
            commandString += command

            if (moveDir): commandString += 'cd %s \n'%(ROOTDIR)
        commandString += bashCheckResultCommand('%s-iteration%d'%(prefix,loop),
                                                jobsPerLoop)
        nJobs -= jobsPerLoop
        idx += jobsPerLoop
        loop += 1

    return commandString

def bashParallelCopyCommand(commands,prefix='job'):
    commandString = ''

    nJobs = len(commands)
    for k in range(nJobs):
        command = ''
        command += '%s &\n' % (commands[k])
        command += 'pids[%d]=$!\n\n' % k
        commandString += command

    commandString += bashCheckResultCommand(prefix,nJobs)

    return commandString

def bashParallelPurgeCommand(filenames,prefix='job'):
    commandString = ''

    nJobs = len(filenames)
    idx = 0
    commandString += 'idx=0\n'
    for k in range(nJobs):
        command = 'if [ -f "%s" ]; then \n' % filenames[k]
        command += 'rm %s &\n' % (filenames[k])
        command += 'pids[${idx}]=$!\n'
        command += 'let "idx+=1"\n'
        command += 'fi\n\n'
        commandString += command
        idx += 1

    commandString += bashCheckResultCommand(prefix,idx)

    return commandString

def purgeDirectoryCommand(baseDirectory):
    bdir = baseDirectory

    commandString = 'shopt -s nullglob\n'
    commandString += 'idx=0\n'

    command = 'for k in {0..%d}\n' % (Nsplit-1)
    command += 'do\n'
    command += '    for file in %s/${k}/*.q\n' % (bdir)
    command += '    do\n'                                                       \
                '       rm ${file} &\n'                                         \
                '       pids[${idx}]=$!\n'                                      \
                '       let "idx+=1"\n'                                         \
                '   done\n'
    command += '    for file in %s/${k}/*.f\n' % (bdir)
    command += '    do\n'                                                       \
                '       rm ${file} &\n'                                         \
                '       pids[${idx}]=$!\n'                                      \
                '       let "idx+=1"\n'                                         \
                '   done\n'
    command += '    for file in %s/${k}/*.dat\n' % (bdir)
    command += '    do\n'                                                       \
                '       rm ${file} &\n'                                         \
                '       pids[${idx}]=$!\n'                                      \
                '       let "idx+=1"\n'                                         \
                '   done\n'
    command += 'done\n'
    commandString += command

    commandString += 'FAIL=0\n'
    commandString += 'k=0\n'
    commandString += 'while [[ $k -lt $idx ]]\n'
    commandString += 'do\n'                                                             \
                     '    wait ${pids[${k}]} || let "FAIL+=1"\n'                        \
                     '    echo "${pids[${k}]}, FAIL: ${FAIL}"\n'                        \
                     '    let "k+=1"\n'                                                 \
                     'done\n\n'
    commandString += 'echo "Number of failures: $FAIL"\n'                               \
                     'if [ $FAIL -ne 0 ]; then\n'
    commandString += '   echo "Purging %s failed."\n' % (bdir)
    commandString += '   exit -1\n'                                                     \
                     'else\n'
    commandString += '   echo "Purging %s succeeded."\n' % (bdir)
    commandString += '   FAIL=0\n'
    commandString += 'fi\n\n'

    # commandString += 'idx=0\n'
    # command = 'for file in %s/* \n' % (DIFFDIR)
    # command += 'do\n'                                                                   \
    #             '   rm ${file} &\n'                                                     \
    #             '   pids[${idx}]=$!\n'                                                  \
    #             '   let "idx+=1"\n'                                                     \
    #             'done\n'
    # commandString += command
    #
    # commandString += 'FAIL=0\n'
    # commandString += 'k=0\n'
    # commandString += 'while [[ $k -lt $idx ]]\n'
    # commandString += 'do\n'                                                             \
    #                  '    wait ${pids[${k}]} || let "FAIL+=1"\n'                        \
    #                  '    echo "${pids[${k}]}, FAIL: ${FAIL}"\n'                        \
    #                  '    let "k+=1"\n'                                                 \
    #                  'done\n\n'
    # commandString += 'echo "Number of failures: $FAIL"\n'                               \
    #                  'if [ $FAIL -ne 0 ]; then\n'
    # commandString += '   echo "Purging diff dir failed."\n'
    # commandString += '   exit -1\n'                                                     \
    #                  'else\n'
    # commandString += '   echo "Purging diff dir succeeded."\n'
    # commandString += '   FAIL=0\n'
    # commandString += 'fi\n\n'

    return commandString
