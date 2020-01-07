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

def QoI(baseDirectory = None):
    if (baseDirectory is None):
        bdir = '.'
    else:
        bdir = baseDirectory

    subJ = np.zeros(3*Nsplit-2)
    J = 0.0
    for k in range(Nsplit):
        subJ[k] = readScalar(bdir + '/' + directories[k] + '/' + outputFiles[k])
        J += subJ[k]
        if (k<Nsplit-1):
            subJ[Nsplit+k] = 0.5 * readScalar(bdir+'/'+diffOutputFiles[k])
            J += matchingConditionWeight * subJ[Nsplit+k]

            subJ[2*Nsplit-1+k] = readScalar(bdir+'/'+lagrangianOutputFiles[k])
            J += subJ[2*Nsplit-1+k]
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

    commands = []
    for j in range(NcontrolRegion):
        commands += ['./zxdoty %s %s %s %s'                                                           \
                    % (outputFiles_[j], xFiles[j], yFiles[j], globalNormFiles[j])]
    commandString += bashParallelLoopCommand(commands,NodesZxdoty,NprocZxdoty,
                                            'zxdoty_control_forcing')

    commands = []
    for k in range(NcontrolRegion,NcontrolRegion+Nsplit-1):
        index = k+1 - NcontrolRegion
        inputFile = '%s/%s'%(directories[index],inputFiles[index])
        commands += ['./spatial_inner_product %s %s --input %s --output %s'                    \
                    % (xFiles[k], yFiles[k], inputFile, outputFiles_[k])]
    commandString += bashParallelLoopCommand(commands,NodesQfileZxdoty,NprocQfileZxdoty,
                                            'zxdoty_initial_condition')

    return commandString

def readInnerProduct(outputFiles_=None):
    if (outputFiles_ is None):
        outputFiles_ = innerProductFiles

    scalar = 0.0
    subScalar = np.zeros(NcontrolRegion+Nsplit-1)
    for k in range(NcontrolRegion+Nsplit-1):
        subScalar[k] = readScalar(outputFiles_[k])
        if (k>=NcontrolRegion):
            ic_idx = k + 1 - NcontrolRegion
            subScalar[k] *= initialConditionControllability[ic_idx]
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

    commands = []
    for j in range(NcontrolRegion):
        command = './zaxpy ' + zFiles[j] + ' ' + "{:.16E}".format(a) + ' ' + xFiles[j]
        if (yFiles is not None):
            command += ' '+yFiles[j]
        commands += [command]
    commandString += bashParallelLoopCommand(commands,NodesZaxpy,NprocZaxpy,
                                            'zaxpy_control_forcing')

    commands = []
    for k in range(NcontrolRegion,NcontrolRegion+Nsplit-1):
        index = k+1 - NcontrolRegion
        w2 = initialConditionControllability[index]
        inputFile = '%s/%s'%(directories[index],inputFiles[index])
        command = './qfile_zaxpy ' + zFiles[k] + ' ' + "{:.16E}".format(a*w2) + ' ' + xFiles[k]
        if (yFiles is not None):
            command += ' '+yFiles[k]
        commands += [command]
    commandString += bashParallelLoopCommand(commands,NodesQfileZaxpy,NprocQfileZaxpy,
                                            'zaxpy_initial_condition')

    return commandString

def distributeCommand(baseDirectory,zeroControlForcing):
    bdir = baseDirectory

    commandString = ''
    commandString += bashGetNodeListCommand()

    if (not zeroControlForcing):
        commands = []
        for k in range(Nsplit):
            kOffset = NtimestepOffset * k
            for j in range(NcontrolRegion):
                sliceControlForcingFile = '%s/%s/%s%s'%(bdir,directories[k],prefixes[k],controlForcingFiles[j])
                globalControlForcingFile = '%s/%s'%(bdir,globalControlSpaceFiles[j])
                commands += ['./slice_control_forcing %s %s %d %d %d'                                                           \
                             % (sliceControlForcingFile, globalControlForcingFile, totalTimestep, kOffset, Nts)]
        commandString += bashParallelLoopCommand(commands,NodesSlice,NprocSlice,'slice_control_forcing')

    for k in range(Nsplit):
        commandString += 'cp ' + '%s/%s'%(bdir,icFiles[k]) + ' ' + '%s/%s/%s'%(bdir,directories[k],icFiles[k]) + '\n'

    return commandString

def gatherControlForcingGradientCommand():
    commandString = ''

    from os import path
    for j in range(NcontrolRegion):
        if (path.exists(globalGradFiles[j])):
            print ("Previous global gradient file %s still exists. Purging it for safety."%(globalGradFiles[j]))
            commandString += 'rm %s \n'%(globalGradFiles[j])

    commands = []
    for k in range(Nsplit):
        kOffset = NtimestepOffset * k
        for j in range(NcontrolRegion):
            sliceGradFile = '%s/%s%s'%(directories[k],prefixes[k],gradientFiles[j])
            commands += ['./paste_control_forcing %s %s %d %d %d'                                                           \
                         % (globalGradFiles[j], sliceGradFile, totalTimestep, kOffset, Nts)]
    commandString += bashParallelLoopCommand(commands,NodesPaste,NprocPaste,'paste_control_forcing')

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

def forwardRunCommand(baseDirectory=None,zeroControlForcing=False):
    if (baseDirectory is None):
        bdir = '.'
    else:
        bdir = baseDirectory

    commandString = ''
    commandString += distributeCommand(bdir,zeroControlForcing)

    commands, commandDirs = [], []
    for k in range(Nsplit):
        commandDirs += ['%s/%s'%(bdir,directories[k])]
        commands += ['./forward --input %s' % inputFiles[k]]
    commandString += bashParallelLoopCommand(commands,NodesForward,NprocForward,
                                            'forward',directories=commandDirs)

    commands = []
    for k in range(Nsplit-1):
        matchingFile = '%s/%s/%s'%(bdir,directories[k],matchingForwardFiles[k])
        icFile = '%s/%s/%s'%(bdir,directories[k+1],icFiles[k+1])
        diffFile = '%s/%s'%(bdir,diffFiles[k])
        inputFile = '%s/%s/%s'%(bdir,directories[k],inputFiles[k])
        commands += ['./qfile_zaxpy %s %.16E %s %s --input %s'                              \
                    % (diffFile,-1.0,icFile,matchingFile,inputFile)]
    commandString += bashParallelLoopCommand(commands,NodesQfileZaxpy,NprocQfileZaxpy,
                                            'forward_qfile_zaxpy')

    commands = []
    for k in range(Nsplit-1):
        diffFile = '%s/%s'%(bdir,diffFiles[k])
        inputFile = '%s/%s/%s'%(bdir,directories[k],inputFiles[k])
        diffOutputFile = '%s/%s'%(bdir,diffOutputFiles[k])
        commands += ['./spatial_inner_product %s %s --output %s --input %s'             \
                    % (diffFile, diffFile, diffOutputFile, inputFile)]
    commandString += bashParallelLoopCommand(commands,NodesQfileZxdoty,NprocQfileZxdoty,
                                            'forward_qfile_zxdoty')

    commands = []
    for k in range(Nsplit-1):
        diffFile = '%s/%s'%(bdir,diffFiles[k])
        inputFile = '%s/%s/%s'%(bdir,directories[k],inputFiles[k])
        lagrangianOutputFile = '%s/%s'%(bdir,lagrangianOutputFiles[k])
        commandString += './spatial_inner_product %s %s --output %s --input %s'             \
                            % (diffFile, lagrangianFiles[k], lagrangianOutputFile, inputFile)
    commandString += bashParallelLoopCommand(commands,NodesQfileZxdoty,NprocQfileZxdoty,
                                            'forward_qfile_lagrangian')

    return commandString

def adjointRunCommand(baseDirectory=None):
    if (baseDirectory is None):
        bdir = '.'
    else:
        bdir = baseDirectory

    commandString = ''
    for k in range(Nsplit-1):
        commandString += 'cd %s/%s \n' % (bdir,directories[k])
        commandString += setOptionCommand(inputFiles[k])
        commandString += 'setOption "enable_adjoint_restart" "true"\n'
        commandString += 'setOption "number_of_timesteps" %d \n' % (Nts-NtimestepOffset)
        commandString += 'setOption "adjoint_restart\/accumulated_timesteps" 0 \n'
        commandString += 'setOption "adjoint_restart\/intermediate_end_timestep" %d \n' % (NtimestepOffset)
        commandString += 'setOption "adjoint_restart\/nonzero_initial_condition" "false" \n'
        if (bdir=='.'):
            commandString += 'cd .. \n'
        else:
            commandString += 'cd ../../ \n'
        commandString += '\n'

    commands, commandDirs = [], []
    for k in range(Nsplit-1):
        commandDirs += ['%s/%s' % (bdir,directories[k])]
        commands += ['./adjoint --input %s' % inputFiles[k]]
    commandString += bashParallelLoopCommand(commands,NodesAdjoint,NprocAdjoint,
                                            'adjoint1',directories=commandDirs)

    k = Nsplit - 1
    commandString += 'cd %s/%s \n' % (bdir,directories[k])
    commandString += setOptionCommand(inputFiles[k])
    commandString += 'setOption "enable_adjoint_restart" "false"\n'
    commandString += 'setOption "number_of_timesteps" %d \n' % (Nts)

    commandString += bashGetNodeListSliceCommand(k,NodesAdjoint)
    commandString += 'srun -N %d -n %d -w ${nodeListString} ' % (NodesAdjoint,NprocAdjoint)
    commandString += './adjoint --input %s' % inputFiles[k]
    commandString += ' &> adjoint_result_%d.out &' % k
    commandString += '\n'
    commandString += 'pids[%d]=$!\n\n' % k
    if (bdir=='.'):
        commandString += 'cd .. \n'
    else:
        commandString += 'cd ../../ \n'
    commandString += '\n'

    commands = []
    for k in range(Nsplit-1):
        matchingAdjointFile = '%s/%s'%(directories[k],matchingAdjointFiles[k])
        inputFile = '%s/%s/%s'%(bdir,directories[k],inputFiles[k])
        diffFile = '%s/%s'%(bdir,diffFiles[k])
        commands += ['./qfile_zaxpy %s %.16E %s %s --input %s'                                      \
        % (matchingAdjointFile, matchingConditionWeight, diffFile, matchingAdjointFile, inputFile)]
    commandString += bashParallelLoopCommand(commands,NodesQfileZaxpy,NprocQfileZaxpy,
                                            'adjoint_run_qfile_zaxpy')

    commands = []
    for k in range(Nsplit-1):
        matchingAdjointFile = '%s/%s'%(directories[k],matchingAdjointFiles[k])
        inputFile = '%s/%s/%s'%(bdir,directories[k],inputFiles[k])
        commands += ['./qfile_zaxpy %s %.16E %s %s --input %s'                                      \
        % (matchingAdjointFile, 1.0, lagrangianFiles[k], matchingAdjointFile, inputFile)]
    commandString += bashParallelLoopCommand(commands,NodesQfileZaxpy,NprocQfileZaxpy,
                                            'adjoint_run_qfile_lagrangian')

    for k in range(Nsplit-1):
        commandString += 'cd %s/%s \n' % (bdir,directories[k])
        commandString += setOptionCommand(inputFiles[k])
        commandString += '\n'
        commandString += 'setOption "number_of_timesteps" %d \n' % (NtimestepOffset)
        commandString += 'setOption "adjoint_restart\/accumulated_timesteps" %d \n' % (Nts-NtimestepOffset)
        commandString += 'setOption "adjoint_restart\/intermediate_end_timestep" 0 \n'
        if (Nts==NtimestepOffset):
            commandString += 'setOption "adjoint_restart\/nonzero_initial_condition" "true" \n'
        if (bdir=='.'):
            commandString += 'cd .. \n'
        else:
            commandString += 'cd ../../ \n'
        commandString += '\n'

    commands, commandDirs = [], []
    for k in range(Nsplit-1):
        commandDirs += ['%s/%s' % (bdir,directories[k])]
        commands += ['./adjoint --input %s' % inputFiles[k]]
    commandString += bashParallelLoopCommand(commands,NodesAdjoint,NprocAdjoint,
                                            'adjoint2',directories=commandDirs)

    # TODO: Either check the last adjoint run, or complete the periodic optimization!!

    for k in range(Nsplit):
        commandString += 'cd %s/%s \n' % (bdir,directories[k])
        commandString += setOptionCommand(inputFiles[k])
        commandString += 'setOption "number_of_timesteps" %d \n' % (Nts)
        commandString += 'setOption "enable_adjoint_restart" "false"\n'
        commandString += 'setOption "adjoint_restart\/nonzero_initial_condition" "false" \n'
        if (bdir=='.'):
            commandString += 'cd .. \n'
        else:
            commandString += 'cd ../../ \n'
        commandString += '\n'

    commands = []
    for k in range(1,Nsplit):
        icAdjointFile = '%s/%s/%s' % (bdir,directories[k],icAdjointFiles[k])
        diffFile = '%s/%s.diff.q'%(bdir,prefixes[k-1])
        icGradFile = '%s/%s' % (bdir,icGradientFiles[k])
        inputFile = '%s/%s/%s'%(bdir,directories[k],inputFiles[k])
        commands += ['./qfile_zaxpy %s %.16E %s %s --input %s'                                      \
                     % (icGradFile, -matchingConditionWeight, diffFile, icAdjointFile, inputFile)]
    commandString += bashParallelLoopCommand(commands,NodesQfileZaxpy,NprocQfileZaxpy,
                                            'adjoint_run_ic')

    commands = []
    for k in range(1,Nsplit):
        icGradFile = '%s/%s' % (bdir,icGradientFiles[k])
        inputFile = '%s/%s/%s'%(bdir,directories[k],inputFiles[k])
        commands += ['./qfile_zaxpy %s %.16E %s %s --input %s'                                     \
                     % (icGradFile, -1.0, lagrangianFiles[k-1], icGradFile, inputFile)]
    commandString += bashParallelLoopCommand(commands,NodesQfileZaxpy,NprocQfileZaxpy,
                                            'adjoint_run_ic_lagrangian')

    return commandString

def dggCommand():
    commandString = ''

    commands = []
    for j in range(NcontrolRegion):
        commands += ['./zwxmwy %s %s %s previous.%s %s'                                                    \
        %(dggFiles[j], globalGradFiles[j], globalGradFiles[j], globalGradFiles[j], globalNormFiles[j])]
    commandString += bashParallelLoopCommand(commands,NodesZxdoty,NprocZxdoty,
                                            'dgg_zwxmwy')

    commands = []
    for k in range(1,Nsplit):
        diffFile = '%s.diff.adjoint.q'%(prefixes[k])
        icGradFile = '%s' % (icGradientFiles[k])
        inputFile = '%s/%s'%(directories[k],inputFiles[k])
        commands += ['./qfile_zaxpy %s %.16E previous.%s %s --input %s'                                    \
                            % (diffFile, -1.0, icGradFile, icGradFile, inputFile)]
    commandString += bashParallelLoopCommand(commands,NodesQfileZaxpy,NprocQfileZaxpy,
                                            'dgg_qfile_zaxpy')

    commands = []
    for k in range(1,Nsplit):
        diffFile = '%s.diff.adjoint.q'%(prefixes[k])
        icGradFile = '%s' % (icGradientFiles[k])
        inputFile = '%s/%s'%(directories[k],inputFiles[k])
        commands += ['./qfile_zaxpy %s %.16E previous.%s %s --input %s'                                    \
                            % (diffFile, -1.0, icGradFile, icGradFile, inputFile)]
    commandString += bashParallelLoopCommand(commands,NodesQfileZaxpy,NprocQfileZaxpy,
                                            'dgg_qfile_zaxpy')

    commands = []
    for k in range(1,Nsplit):
        idx = NcontrolRegion-1 + k
        diffFile = '%s.diff.adjoint.q'%(prefixes[k])
        inputFile = '%s/%s'%(directories[k],inputFiles[k])
        commands += ['./spatial_inner_product %s %s --output %s --input %s'                                \
                            % (diffFile, icGradientFiles[k], dggFiles[idx], inputFile)]
    commandString += bashParallelLoopCommand(commands,NodesQfileZxdoty,NprocQfileZxdoty,
                                            'dgg_qfile_zxdoty')

    return commandString
