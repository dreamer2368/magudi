from constants import *
from filenames import *

import numpy as np
import subprocess
import pandas as pd

__all__ = ['commandScriptor','bashCommandScriptor','fluxCommandScriptor','scriptorSwitcher']

class commandScriptor:
    def singleJobCommand(self,commands,procedure,prefix='job',directories=None):
        pass

    def parallelLoopCommand(self,commands,procedure,prefix='job',directories=None):
        pass

    def nonMPILoopCommand(self,commands,prefix='job'):
        commandString = ''

        nJobs = len(commands)
        for k in range(nJobs):
            command = ''
            command += '%s &\n' % (commands[k])
            command += 'pids[%d]=$!\n\n' % k
            commandString += command

        commandString += self.checkResultCommand(prefix,nJobs)

        return commandString

    def parallelPurgeCommand(self,filenames,prefix='job'):
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

        commandString += self.checkResultCommand(prefix,idx)

        return commandString

    def purgeDirectoryCommand(self,baseDirectory):
        bdir = baseDirectory

        commandString = 'shopt -s nullglob\n'
        commandString += 'idx=0\n'

        command = 'for file in %s/**/*.q\n' % (bdir)
        command += 'do\n'                                                       \
                   '   rm ${file} &\n'                                         \
                   '   pids[${idx}]=$!\n'                                      \
                   '   let "idx+=1"\n'                                         \
                   'done\n'
        command += 'for file in %s/**/*.f\n' % (bdir)
        command += 'do\n'                                                       \
                   '    rm ${file} &\n'                                         \
                   '    pids[${idx}]=$!\n'                                      \
                   '    let "idx+=1"\n'                                         \
                   'done\n'
        command += 'for file in %s/**/*.dat\n' % (bdir)
        command += 'do\n'                                                       \
                   '    rm ${file} &\n'                                         \
                   '    pids[${idx}]=$!\n'                                      \
                   '    let "idx+=1"\n'                                         \
                   'done\n'
        commandString += command

        commandString += 'FAIL=0\n'
        commandString += 'k=0\n'
        commandString += 'while [[ $k -lt $idx ]]\n'
        commandString += 'do\n'                                                             \
                         '    wait ${pids[${k}]} || let "FAIL+=1"\n'                        \
                         '    echo "${k}, FAIL: ${FAIL}"\n'                                 \
                         '    let "k+=1"\n'                                                 \
                         'done\n\n'
        commandString += 'echo "Number of failures: $FAIL"\n'                               \
                         'shopt -u nullglob\n'                                              \
                         'if [ $FAIL -eq 0 ]; then\n'
        commandString += '   echo "Purging %s succeeded."\n' % (bdir)
        commandString += '   FAIL=0\n'                                                      \
                         'else\n'
        commandString += '   numFiles=0\n'                                                  \
                         '   let "numFiles+=$(ls %s/**/*.q | wc -l)"\n' % (bdir)
        commandString += '   let "numFiles+=$(ls %s/**/*.f | wc -l)"\n' % (bdir)
        commandString += '   let "numFiles+=$(ls %s/**/*.dat | wc -l)"\n' % (bdir)
        commandString += '   if [ $numFiles -eq 0 ]; then\n'
        commandString += '      echo "No file remains in subdirectories. Purging %s succeeded."\n' % (bdir)
        commandString += '      FAIL=0\n'                                                   \
                         '   else\n'                                                        \
                         '      echo "Purging %s failed."\n' % (bdir)
        commandString += '      exit -1\n'                                                  \
                         '   fi\n'
        commandString += 'fi\n\n'

        return commandString

    def checkResultCommand(self,procedureName,countFails=0):
        if(countFails>0):
            commandString = 'FAIL=0\n'
            commandString += 'for k in {0..%d}\n'%(countFails-1)
            commandString += 'do\n'                                                             \
                             '    wait ${pids[${k}]} || let "FAIL+=1"\n'                        \
                             '    echo "${k}, FAIL: ${FAIL}"\n'                                 \
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

class bashCommandScriptor(commandScriptor):
    def getNodeListCommand(self):
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

    def getNodeListSliceCommand(self,index, numNodes):
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

    def singleJobCommand(self,commands,procedure,prefix='job',directories=None):
        if (directories is None):
            moveDir = False
        else:
            moveDir = True
            if (len(commands)!=len(directories)):
                raise ValueError('Provide directories for all commands.')

        nodePerCommand, procsPerCommand = procedureSwitcher.get(procedure)

        commandString = ''
        commandString += self.getNodeListCommand()

        for k, command in enumerate(commands):
            if (moveDir): commandString += 'cd %s \n'%(directories[k])
            commandString += self.getNodeListSliceCommand(0,nodePerCommand)
            commandString += 'srun -N %d -n %d -w ${nodeListString} '                         \
                             % (nodePerCommand,procsPerCommand)
            commandString += command
            commandString += ' &> %s/%s_result_%d.out \n' % (OUTDIR,prefix,k)
            commandString += self.checkResultCommand('%s-serial-%d'%(prefix,k))
            if (moveDir): commandString += 'cd %s \n'%(ROOTDIR)

        return commandString

    def parallelLoopCommand(self,commands,procedure,prefix='job',directories=None):
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

        if (enableParallelBash): commandString += self.getNodeListCommand()
        idx, nJobs = 0, len(commands)
        loop = 0
        while (nJobs>0):
            jobsPerLoop = maxJobsPerLoop if (nJobs>maxJobsPerLoop) else nJobs
            for k in range(jobsPerLoop):
                if (moveDir): commandString += 'cd %s \n'%(directories[idx+k])

                command = ''
                if (enableParallelBash):
                    command += self.getNodeListSliceCommand(k,nodePerCommand)
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
            commandString += self.checkResultCommand('%s-iteration%d'%(prefix,loop),
                                                    jobsPerLoop)
            nJobs -= jobsPerLoop
            idx += jobsPerLoop
            loop += 1

        return commandString

class fluxCommandScriptor(commandScriptor):
    def waitallCommand(self):
        commandString = 'waitall() {\n'                                                     \
                        'local rc=0\n'                                                      \
                        'local rcsum=0\n'                                                   \
                        'flux queue drain\n'                                                \
                        'local k=0\n'                                                       \
                        'for job in $@\n'                                                   \
                        'do\n'                                                              \
                        '    flux job status ${job} || rc=$?\n'                             \
                        '    echo "${k}-th job id: ${job}, FAIL: ${rc}"\n'                  \
                        '    if [ $rc -ne 0 ]; then\n'                                      \
                        '        flux job attach ${job}\n'                                  \
                        '        flux job info ${job} R\n'                                  \
                        '        let "rcsum+=1"\n'                                          \
                        '    fi\n'                                                          \
                        '    let "k+=1"; rc=0\n'                                            \
                        'done\n'                                                            \
                        'if [ $rcsum -gt 0 ]; then\n'                                       \
                        '    return 1\n'                                                    \
                        'else\n'                                                            \
                        '    return 0\n'                                                    \
                        'fi\n'                                                              \
                        '}\n\n'
        return commandString

    def singleJobCommand(self,commands,procedure,prefix='job',directories=None):
        if (directories is None):
            moveDir = False
        else:
            moveDir = True
            if (len(commands)!=len(directories)):
                raise ValueError('Provide directories for all commands.')

        nodePerCommand, procsPerCommand = procedureSwitcher.get(procedure)

        commandString = self.waitallCommand()
        commandString += 'JOBIDS=""\n'

        for k, command in enumerate(commands):
            if (moveDir): commandString += 'cd %s \n'%(directories[k])
            commandString += 'JOBIDS="$JOBIDS $(flux mini submit -n %d --output=%s/%s_result_%d.out '                         \
                             % (procsPerCommand,OUTDIR,prefix,k)
            commandString += command
            commandString += ')"\n\nwaitall ${JOBIDS}\n'
            commandString += self.checkResultCommand('%s-serial-%d'%(prefix,k))
            if (moveDir): commandString += 'cd %s \n'%(ROOTDIR)

        return commandString

    def parallelLoopCommand(self,commands,procedure,prefix='job',directories=None):
        if (directories is None):
            moveDir = False
        else:
            moveDir = True
            if (len(commands)!=len(directories)):
                raise ValueError('Provide directories for all commands.')

        nodePerCommand, procsPerCommand = procedureSwitcher.get(procedure)

        maxJobsPerLoop = int(np.floor(maxProcs/procsPerCommand))
        if (maxJobsPerLoop<1):
            raise ValueError('The job %s requires at least %d number of procs. '                \
                                %(prefix,procsPerCommand) +                                      \
                             'Only %d procs are assigned.'%(maxProcs))

        commandString = self.waitallCommand()

        idx, nJobs = 0, len(commands)
        loop = 0
        while (nJobs>0):
            commandString += 'JOBIDS=""\n'
            jobsPerLoop = maxJobsPerLoop if (nJobs>maxJobsPerLoop) else nJobs
            for k in range(jobsPerLoop):
                if (moveDir): commandString += 'cd %s \n'%(directories[idx+k])

                command = ''
                command += 'JOBIDS="$JOBIDS $(flux mini submit -n %d --output=%s/%s_result_%d.out '                         \
                             % (procsPerCommand,OUTDIR,prefix,idx+k)
                command += commands[idx+k]
                command += ')"\n'
                commandString += command

                if (moveDir): commandString += 'cd %s \n'%(ROOTDIR)
            commandString += 'waitall ${JOBIDS}\n'
            commandString += self.checkResultCommand('%s-iteration%d'%(prefix,loop))
            nJobs -= jobsPerLoop
            idx += jobsPerLoop
            loop += 1

        return commandString

scriptorSwitcher = {
    'bash': bashCommandScriptor(),
    'flux': fluxCommandScriptor()
}
