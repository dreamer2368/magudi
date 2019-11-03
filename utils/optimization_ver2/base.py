from constants import *
from filenames import *
from parallelBashCommand import *

import numpy as np
import subprocess
import pandas as pd

magudiSetOptionCommand = 'function setOption() {\n'                                         \
                         '    if grep -q "$1" magudi.inp\n'                                 \
                         '    then\n'                                                       \
                         '	sed -i "s/^.*$1.*$/$1 = $2/g" magudi.inp\n'                     \
                         '    else\n'                                                       \
                         '	echo "$1 = $2" >> magudi.inp\n'                                 \
                         '    fi\n'                                                         \
                         '}\n'

def setOptionCommand(inputFilename):
    commandString = 'function setOption() {\n'                                          \
                     '    if grep -q "$1" ' + inputFilename + '\n'                      \
                     '    then\n'                                                       \
                     '    sed -i.bu "s/^.*$1.*$/$1 = $2/g" ' + inputFilename + '\n'     \
                     '    else\n'                                                       \
                     '    echo "$1 = $2" >> ' + inputFilename + '\n'                    \
                     '    fi\n'                                                         \
                     '}\n'
    return commandString

def readScalar(scalarFilename):
    try:
        fID = open(scalarFilename,'r')
    except FileNotFoundError:
        print (scalarFilename+' is not found.')
        return np.nan
    else:
        scalar = float(fID.read())
        fID.close()
        return scalar

def QoI(rootDirectory = None):
    if (rootDirectory is None):
        rdir = '.'
    else:
        rdir = rootDirectory

    subJ = np.zeros(2*Nsplit-1)
    J = 0.0
    for k in range(Nsplit):
        subJ[k] = readScalar(rdir + '/' + directories[k] + '/' + outputFiles[k])
        J += subJ[k]
        if (k<Nsplit-1):
            subJ[Nsplit+k] = 0.5 * readScalar(rdir+'/'+diffOutputFiles[k])
            J += matchingConditionWeight * subJ[Nsplit+k]
    return J, subJ

def collectQoI(logFilename):
    import pandas as pd

    new_J, dummy = QoI('x')

    df = pd.read_csv(logFilename, sep='\t', header=0)
    df.at[df['directory index']=='x','QoI'] = new_J
    df = df.sort_values(by='step',ascending=True)
    df = df.reset_index(drop=True)
    return df

def innerProductCommand(xFiles,yFiles,outputFiles_=None):
    if (outputFiles_ is None):
        outputFiles_ = innerProductFiles

    Nx, Ny = len(xFiles), len(yFiles)
    if ((Nx!=Nsplit+NcontrolRegion-1) or (Ny!=Nsplit+NcontrolRegion-1)):
        raise LookupError('Number of files must be equal to '+str(Nsplit+NcontrolRegion-1)+'!')

    commandString = ''
    commandString += bashGetNodeListCommand()
    for j in range(NcontrolRegion):
        command = bashGetNodeListSliceCommand(j,NodesZxdoty)
        command += 'srun -N %d -n %d -w ${nodeListString} ' % (NodesZxdoty,NprocZxdoty)
        command += './zxdoty %s %s %s %s'                                                           \
                    % (outputFiles_[j], xFiles[j], yFiles[j], globalNormFiles[j])
        command += ' &> zxdoty_result_%d.out &' % j
        command += '\n'
        command += 'pids[%d]=$!\n\n' % j
        commandString += command

    commandString += bashCheckResultCommand('zxdoty_control_forcing',NcontrolRegion)

    for k in range(NcontrolRegion,NcontrolRegion+Nsplit-1):
        index = k+1 - NcontrolRegion
        command = bashGetNodeListSliceCommand(index,NodesQfileZxdoty)
        command += 'srun -N %d -n %d -w ${nodeListString} ' % (NodesQfileZxdoty,NprocQfileZxdoty)

        inputFile = '%s/%s'%(directories[index],inputFiles[index])
        command += './spatial_inner_product %s %s --input %s --output %s'                    \
                    % (xFiles[k], yFiles[k], inputFile, outputFiles_[k])
        command += ' &> zxdoty_result_%d.out &' % k
        command += '\n'
        command += 'pids[%d]=$!\n\n' % k
        commandString += command

    commandString += bashCheckResultCommand('zxdoty_initial_condition',Nsplit-1)

    return commandString

def readInnerProduct(outputFiles_=None):
    if (outputFiles_ is None):
        outputFiles_ = innerProductFiles

    scalar = 0.0
    subScalar = np.zeros(NcontrolRegion+Nsplit-1)
    for k in range(NcontrolRegion+Nsplit-1):
        subScalar[k] = readScalar(outputFiles_[k])
        scalar += subScalar[k]

    return scalar, subScalar

def zaxpyCommand(zFiles,a,xFiles,yFiles=None):
    Nx, Nz = len(xFiles), len(zFiles)
    if ((Nx!=Nsplit+NcontrolRegion-1) or (Nz!=Nsplit+NcontrolRegion-1)):
        raise LookupError('Number of files must be equal to '+str(Nsplit+NcontrolRegion-1)+'!')
    if (yFiles is not None):
        Ny = len(yFiles)
        if (Ny!=Nsplit+NcontrolRegion-1):
            raise LookupError('Number of files must be equal to '+str(Nsplit+NcontrolRegion-1)+'!')

    commandString = ''
    commandString += bashGetNodeListCommand()
    for j in range(NcontrolRegion):
        commandString += bashGetNodeListSliceCommand(j,NodesZaxpy)
        commandString += 'srun -N %d -n %d -w ${nodeListString} ' % (NodesZaxpy,NprocZaxpy)

        commandString += './zaxpy ' + zFiles[j] + ' ' + "{:.16E}".format(a) + ' ' + xFiles[j]
        if (yFiles is not None):
            commandString += ' '+yFiles[j]
        commandString += ' &> zaxpy_result_%d.out &' % j
        commandString += '\n'
        commandString += 'pids[%d]=$!\n\n' % j

    commandString += bashCheckResultCommand('zaxpy_control_forcing',NcontrolRegion)

    w2 = initialConditionControllability
    for k in range(NcontrolRegion,NcontrolRegion+Nsplit-1):
        index = k+1 - NcontrolRegion
        commandString += bashGetNodeListSliceCommand(index,NodesQfileZaxpy)
        commandString += 'srun -N %d -n %d -w ${nodeListString} ' % (NodesQfileZaxpy,NprocQfileZaxpy)

        inputFile = '%s/%s'%(directories[index],inputFiles[index])
        commandString += './qfile_zaxpy ' + zFiles[k] + ' ' + "{:.16E}".format(a*w2) + ' ' + xFiles[k]
        if (yFiles is not None):
            commandString += ' '+yFiles[k]
        commandString += ' &> zaxpy_result_%d.out &' % k
        commandString += '\n'
        commandString += 'pids[%d]=$!\n\n' % k

    commandString += bashCheckResultCommand('zaxpy_initial_condition',Nsplit-1)

    return commandString

def distributeCommand(rootDirectory):
    rdir = rootDirectory

    commandString = ''
    commandString += bashGetNodeListCommand()
    for k in range(Nsplit):
        kOffset = NtimestepOffset * k
        for j in range(NcontrolRegion):
            nodeIndex = NcontrolRegion * k + j
            commandString += bashGetNodeListSliceCommand(nodeIndex,NodesSlice)
            commandString += 'srun -N %d -n %d -w ${nodeListString} ' % (NodesSlice,NprocSlice)

            sliceControlForcingFile = '%s/%s/%s%s'%(rdir,directories[k],prefixes[k],controlForcingFiles[j])
            globalControlForcingFile = '%s/%s'%(rdir,globalControlSpaceFiles[j])
            commandString += './slice_control_forcing %s %s %d %d %d'                                                           \
                             % (sliceControlForcingFile, globalControlForcingFile, totalTimestep, kOffset, Nts)
            commandString += ' &> slice_result_%d.out &' % nodeIndex
            commandString += '\n'
            commandString += 'pids[%d]=$!\n\n' % nodeIndex

    commandString += bashCheckResultCommand('slice_control_forcing',Nsplit*NcontrolRegion)

    for k in range(Nsplit):
        commandString += 'cp ' + '%s/%s'%(rdir,icFiles[k]) + ' ' + '%s/%s/%s'%(rdir,directories[k],icFiles[k]) + '\n'

    return commandString

def gatherControlForcingGradientCommand():
    commandString = ''
    commandString += bashGetNodeListCommand()
    for k in range(Nsplit):
        kOffset = NtimestepOffset * k
        for j in range(NcontrolRegion):
            commandString += bashGetNodeListSliceCommand(j,NodesPaste)
            commandString += 'srun -N %d -n %d -w ${nodeListString} ' % (NodesPaste,NprocPaste)

            sliceGradFile = '%s/%s%s'%(directories[k],prefixes[k],gradientFiles[j])
            commandString += './paste_control_forcing %s %s %d %d %d'                                                           \
                             % (globalGradFiles[j], sliceGradFile, totalTimestep, kOffset, Nts)
            commandString += ' &> paste_result_%d.out &' % j
            commandString += '\n'
            commandString += 'pids[%d]=$!\n\n' % j
        commandString += bashCheckResultCommand('paste_control_forcing_%d'%k,NcontrolRegion)
    return commandString

def switchDirectory(firstDirectory, secondDirectory, df=None):
    if (firstDirectory==secondDirectory):
        return

    import subprocess
    subprocess.check_call('mv '+str(firstDirectory)+' temp', shell=True)
    subprocess.check_call('mv '+str(secondDirectory)+' '+str(firstDirectory), shell=True)
    subprocess.check_call('mv temp '+str(secondDirectory), shell=True)

    import pandas as pd
    if isinstance(df,pd.DataFrame):
        df.at[df['directory index']==firstDirectory,'directory index'] = -1
        df.at[df['directory index']==secondDirectory,'directory index'] = firstDirectory
        df.at[df['directory index']==-1,'directory index'] = secondDirectory
    return

def forwardRunCommand(rootDirectory=None):
    if (rootDirectory is None):
        rdir = '.'
    else:
        rdir = rootDirectory

    commandString = ''
    commandString += distributeCommand(rdir)

    for k in range(Nsplit):
        commandString += 'cd %s/%s \n'%(rdir,directories[k])

        commandString += bashGetNodeListSliceCommand(k,NodesForward)
        commandString += 'srun -N %d -n %d -w ${nodeListString} ' % (NodesForward,NprocForward)
        commandString += './forward --input %s' % inputFiles[k]
        commandString += ' &> forward_result_%d.out &' % k
        commandString += '\n'
        commandString += 'pids[%d]=$!\n\n' % k

        if (rdir=='.'):
            commandString += 'cd .. \n'
        else:
            commandString += 'cd ../../ \n'

    commandString += bashCheckResultCommand('forward_runs',Nsplit)

    for k in range(Nsplit-1):
        commandString += bashGetNodeListSliceCommand(k,NodesQfileZaxpy)
        commandString += 'srun -N %d -n %d -w ${nodeListString} ' % (NodesQfileZaxpy,NprocQfileZaxpy)

        kOffset = startTimestep + NtimestepOffset * (k+1)
        matchingFile = '%s/%s/%s'%(rdir,directories[k],matchingForwardFiles[k])
        icFile = '%s/%s/%s'%(rdir,directories[k+1],icFiles[k+1])
        diffFile = '%s/%s'%(rdir,diffFiles[k])
        inputFile = '%s/%s/%s'%(rdir,directories[k],inputFiles[k])
        commandString += './qfile_zaxpy %s %.16E %s %s --input %s'                          \
                            % (diffFile,-1.0,icFile,matchingFile,inputFile)
        commandString += ' &> forward_qfile_zaxpy_result_%d.out &' % k
        commandString += '\n'
        commandString += 'pids[%d]=$!\n\n' % k

    commandString += bashCheckResultCommand('forward_run_qfile_zaxpy',Nsplit-1)

    for k in range(Nsplit-1):
        commandString += bashGetNodeListSliceCommand(k,NodesQfileZxdoty)
        commandString += 'srun -N %d -n %d -w ${nodeListString} ' % (NodesQfileZxdoty,NprocQfileZxdoty)

        diffFile = '%s/%s'%(rdir,diffFiles[k])
        inputFile = '%s/%s/%s'%(rdir,directories[k],inputFiles[k])
        diffOutputFile = '%s/%s'%(rdir,diffOutputFiles[k])
        commandString += './spatial_inner_product %s %s --output %s --input %s'             \
                            % (diffFile, diffFile, diffOutputFile, inputFile)
        commandString += ' &> forward_qfile_zxdoty_result_%d.out &' % k
        commandString += '\n'
        commandString += 'pids[%d]=$!\n\n' % k

    commandString += bashCheckResultCommand('forward_run_qfile_zxdoty',Nsplit-1)

    return commandString

def adjointRunCommand(rootDirectory=None):
    if (rootDirectory is None):
        rdir = '.'
    else:
        rdir = rootDirectory

    commandString = ''
    commandString += bashGetNodeListCommand()
    for k in range(Nsplit-1):
        commandString += 'cd %s/%s \n' % (rdir,directories[k])
        commandString += setOptionCommand(inputFiles[k])
        commandString += 'setOption "enable_adjoint_restart" "true"\n'
        commandString += 'setOption "number_of_timesteps" %d \n' % (Nts-NtimestepOffset)
        commandString += 'setOption "adjoint_restart\/accumulated_timesteps" 0 \n'
        commandString += 'setOption "adjoint_restart\/intermediate_end_timestep" %d \n' % (NtimestepOffset)

        commandString += bashGetNodeListSliceCommand(k,NodesAdjoint)
        commandString += 'srun -N %d -n %d -w ${nodeListString} ' % (NodesAdjoint,NprocAdjoint)
        commandString += './adjoint --input %s' % inputFiles[k]
        commandString += ' &> adjoint_result_%d.out &' % k
        commandString += '\n'
        commandString += 'pids[%d]=$!\n\n' % k

        if (rdir=='.'):
            commandString += 'cd .. \n'
        else:
            commandString += 'cd ../../ \n'
        commandString += '\n'

    k = Nsplit - 1
    commandString += 'cd %s/%s \n' % (rdir,directories[k])
    commandString += setOptionCommand(inputFiles[k])
    commandString += 'setOption "enable_adjoint_restart" "false"\n'
    commandString += 'setOption "number_of_timesteps" %d \n' % (Nts)

    commandString += bashGetNodeListSliceCommand(k,NodesAdjoint)
    commandString += 'srun -N %d -n %d -w ${nodeListString} ' % (NodesAdjoint,NprocAdjoint)
    commandString += './adjoint --input %s' % inputFiles[k]
    commandString += ' &> adjoint_result_%d.out &' % k
    commandString += '\n'
    commandString += 'pids[%d]=$!\n\n' % k
    if (rdir=='.'):
        commandString += 'cd .. \n'
    else:
        commandString += 'cd ../../ \n'
    commandString += '\n'

    commandString += bashCheckResultCommand('adjoint_run_1',Nsplit-1)

    for k in range(Nsplit-1):
        commandString += bashGetNodeListSliceCommand(k,NodesQfileZaxpy)
        commandString += 'srun -N %d -n %d -w ${nodeListString} ' % (NodesQfileZaxpy,NprocQfileZaxpy)

        matchingAdjointFile = '%s/%s'%(directories[k],matchingAdjointFiles[k])
        inputFile = '%s/%s/%s'%(rdir,directories[k],inputFiles[k])
        diffFile = '%s/%s'%(rdir,diffFiles[k])
        commandString += './qfile_zaxpy %s %.16E %s %s --input %s'                                      \
        % (matchingAdjointFile, matchingConditionWeight, diffFile, matchingAdjointFile, inputFile)
        commandString += ' &> adjoint_qfile_zaxpy_result_%d.out &' % k
        commandString += '\n'
        commandString += 'pids[%d]=$!\n\n' % k

    commandString += bashCheckResultCommand('adjoint_run_qfile_zaxpy',Nsplit-1)

    for k in range(Nsplit-1):
        commandString += 'cd %s/%s \n' % (rdir,directories[k])
        commandString += setOptionCommand(inputFiles[k])
        commandString += '\n'
        commandString += 'setOption "number_of_timesteps" %d \n' % (NtimestepOffset)
        commandString += 'setOption "adjoint_restart\/accumulated_timesteps" %d \n' % (Nts-NtimestepOffset)
        commandString += 'setOption "adjoint_restart\/intermediate_end_timestep" 0 \n'

        commandString += bashGetNodeListSliceCommand(k,NodesAdjoint)
        commandString += 'srun -N %d -n %d -w ${nodeListString} ' % (NodesAdjoint,NprocAdjoint)
        commandString += './adjoint --input %s' % inputFiles[k]
        commandString += ' &> adjoint_result_%d.out &' % k
        commandString += '\n'
        commandString += 'pids[%d]=$!\n\n' % k

        commandString += 'setOption "number_of_timesteps" %d \n' % (Nts)
        commandString += 'setOption "enable_adjoint_restart" "false"\n'
        if (rdir=='.'):
            commandString += 'cd .. \n'
        else:
            commandString += 'cd ../../ \n'
        commandString += '\n'

    commandString += bashCheckResultCommand('adjoint_run_2',Nsplit)

    for k in range(1,Nsplit):
        commandString += bashGetNodeListSliceCommand(k-1,NodesQfileZaxpy)
        commandString += 'srun -N %d -n %d -w ${nodeListString} ' % (NodesQfileZaxpy,NprocQfileZaxpy)

        icAdjointFile = '%s/%s/%s' % (rdir,directories[k],icAdjointFiles[k])
        diffFile = '%s/%s.diff.q'%(rdir,prefixes[k-1])
        icGradFile = '%s/%s' % (rdir,icGradientFiles[k])
        inputFile = '%s/%s/%s'%(rdir,directories[k],inputFiles[k])
        commandString += './qfile_zaxpy %s %.16E %s %s --input %s'                                \
                         % (icGradFile, -matchingConditionWeight, diffFile, icAdjointFile, inputFile)
        commandString += ' &> adjoint_initial_condition_%d.out &' % (k-1)
        commandString += '\n'
        commandString += 'pids[%d]=$!\n\n' % (k-1)

    commandString += bashCheckResultCommand('adjoint_run_ic',Nsplit-1)

    return commandString

def dggCommand():
    commandString = ''
    commandString += bashGetNodeListCommand()
    for j in range(NcontrolRegion):
        commandString += bashGetNodeListSliceCommand(j,NodesZxdoty)
        commandString += 'srun -N %d -n %d -w ${nodeListString} ' % (NodesZxdoty,NprocZxdoty)
        commandString += './zwxmwy %s %s %s previous.%s %s'                                               \
        %(dggFiles[j], globalGradFiles[j], globalGradFiles[j], globalGradFiles[j], globalNormFiles[j])
        commandString += ' &> dgg_zwxmwy_%d.out &' % j
        commandString += '\n'
        commandString += 'pids[%d]=$!\n\n' % j

    commandString += bashCheckResultCommand('dgg_zwxmwy',NcontrolRegion)

    for k in range(1,Nsplit):
        commandString += bashGetNodeListSliceCommand(k-1,NodesQfileZaxpy)
        commandString += 'srun -N %d -n %d -w ${nodeListString} ' % (NodesQfileZaxpy,NprocQfileZaxpy)

        diffFile = '%s.diff.adjoint.q'%(prefixes[k])
        icGradFile = '%s' % (icGradientFiles[k])
        inputFile = '%s/%s'%(directories[k],inputFiles[k])
        commandString += './qfile_zaxpy %s %.16E previous.%s %s --input %s'                                \
                            % (diffFile, -1.0, icGradFile, icGradFile, inputFile)
        commandString += ' &> dgg_qfile_zaxpy_%d.out &' % (k-1)
        commandString += '\n'
        commandString += 'pids[%d]=$!\n\n' % (k-1)

    commandString += bashCheckResultCommand('dgg_qfile_zaxpy',Nsplit-1)

    for k in range(1,Nsplit):
        commandString += bashGetNodeListSliceCommand(k-1,NodesQfileZxdoty)
        commandString += 'srun -N %d -n %d -w ${nodeListString} ' % (NodesQfileZxdoty,NprocQfileZxdoty)

        idx = NcontrolRegion-1 + k
        diffFile = '%s.diff.adjoint.q'%(prefixes[k])
        inputFile = '%s/%s'%(directories[k],inputFiles[k])
        commandString += './spatial_inner_product %s %s --output %s --input %s'                            \
                            % (diffFile, icGradientFiles[k], dggFiles[idx], inputFile)
        commandString += ' &> dgg_qfile_zxdoty_%d.out &' % (k-1)
        commandString += '\n'
        commandString += 'pids[%d]=$!\n\n' % (k-1)

    commandString += bashCheckResultCommand('dgg_qfile_zxdoty',Nsplit-1)

    return commandString
