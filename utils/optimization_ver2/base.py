from constants import *
from filenames import *
from command_scriptor import scriptorSwitcher

import numpy as np
import subprocess
import pandas as pd

scriptor = scriptorSwitcher.get(scriptorType)

magudiSetOptionCommand = 'function setOption() {\n'                                         \
                         '    if grep -q "$1" magudi.inp\n'                                 \
                         '    then\n'                                                       \
                         '	sed -i "s+^.*$1.*$+$1 = $2+g" magudi.inp\n'                     \
                         '    else\n'                                                       \
                         '	echo "$1 = $2" >> magudi.inp\n'                                 \
                         '    fi\n'                                                         \
                         '}\n'

generalSetOptionCommand = 'cat <<EOF > setOption.sh\n'                                   \
                          '#!/bin/bash\n'                                                \
                          'if grep -q "\$2" \$1 \n'                                      \
                          'then\n'                                                       \
                          'sed -i.bu "s+^.*\$2.*\$+\$2 = \$3+g" \$1 \n'                  \
                          'else\n'                                                       \
                          'echo "\$2 = \$3" >> \$1 \n'                                   \
                          'fi\n'                                                         \
                          'EOF\n'                                                        \
                          'chmod u+x setOption.sh\n'

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

def QoI(baseDirectory = 'x0'):
    bdir = baseDirectory

    subJ = np.zeros(3*Nsplit)
    J = 0.0
    for k in range(Nsplit):
        subJ[k] = readScalar(bdir + '/' + directories[k] + '/' + outputFiles[k])
        if (not ignoreObjective):
            J += subJ[k]

        subJ[Nsplit+k] = 0.5 * readScalar(diffOutputFiles[k])
        J += matchingConditionWeight[k] * subJ[Nsplit+k]

        if (useLagrangian):
            subJ[2*Nsplit+k] = readScalar(lagrangianOutputFiles[k])
            # followed the convention of Portryagin's minimum principle
            J -= subJ[2*Nsplit+k]
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
    if ((Nx!=NcontrolSpace) or (Ny!=NcontrolSpace)):
        raise LookupError('Number of files must be equal to '+str(NcontrolSpace)+'!')

    commandString = ''

    commands = []
    for j in range(NcontrolRegion):
        commands += ['./zxdoty %s %s %s %s'                                                           \
                    % (outputFiles_[j], xFiles[j], yFiles[j], globalNormFiles[j])]
    commandString += scriptor.parallelLoopCommand(commands,'zxdoty','zxdoty_control_forcing')

    commands = []
    for k in range(NcontrolRegion,NcontrolSpace):
        index = k - NcontrolRegion
        commands += ['./spatial_inner_product %s %s --input %s --output %s'                    \
                    % (xFiles[k], yFiles[k], globalInputFile, outputFiles_[k])]
    commandString += scriptor.parallelLoopCommand(commands,'qfile-zxdoty',
                                                  'zxdoty_initial_condition')

    return commandString

def readInnerProduct(outputFiles_=None):
    if (outputFiles_ is None):
        outputFiles_ = innerProductFiles

    scalar = 0.0
    subScalar = np.zeros(NcontrolSpace)
    for k in range(NcontrolSpace):
        subScalar[k] = readScalar(outputFiles_[k])
    subScalar[-Nsplit:] *= initialConditionControllability
    scalar = np.sum(subScalar)

    return scalar, subScalar

def zaxpyCommand(zFiles,a,xFiles,yFiles=None):
    Nx, Nz = len(xFiles), len(zFiles)
    if ((Nx!=NcontrolSpace) or (Nz!=NcontrolSpace)):
        raise LookupError('Number of files must be equal to '+str(NcontrolSpace)+'!')
    if (yFiles is not None):
        Ny = len(yFiles)
        if (Ny!=NcontrolSpace):
            raise LookupError('Number of files must be equal to '+str(NcontrolSpace)+'!')

    commandString = ''

    commands = []
    for j in range(NcontrolRegion):
        command = './zaxpy ' + zFiles[j] + ' ' + "{:.16E}".format(a) + ' ' + xFiles[j]
        if (yFiles is not None):
            command += ' '+yFiles[j]
        commands += [command]
    commandString += scriptor.parallelLoopCommand(commands,'zaxpy',
                                                  'zaxpy_control_forcing')

    commands = []
    for k in range(NcontrolRegion,NcontrolSpace):
        index = k - NcontrolRegion
        w2 = initialConditionControllability[index]
        command = './qfile_zaxpy ' + zFiles[k] + ' ' + "{:.16E}".format(a*w2) + ' ' + xFiles[k]
        if (yFiles is not None):
            command += ' '+yFiles[k]
        commands += [command]
    commandString += scriptor.parallelLoopCommand(commands,'qfile-zaxpy',
                                                  'zaxpy_initial_condition')

    return commandString

def distributeCommand(baseDirectory,zeroControlForcing):
    bdir = baseDirectory

    commandString = ''

    if (not zeroControlForcing):
        commands = []
        for k in range(Nsplit):
            kOffset = Nts * k
            for j in range(NcontrolRegion):
                sliceControlForcingFile = '%s/%s/%s%s'%(bdir,directories[k],prefixes[k],controlForcingFiles[j])
                globalControlForcingFile = '%s/%s'%(bdir,globalControlSpaceFiles[j])
                commands += ['./slice_control_forcing %s %s %d %d %d'                                                           \
                             % (sliceControlForcingFile, globalControlForcingFile, totalTimestep, kOffset, Nts)]
        commandString += scriptor.parallelLoopCommand(commands,'slice','slice_control_forcing')

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
        kOffset = Nts * k
        for j in range(NcontrolRegion):
            sliceGradFile = 'x0/%s/%s%s'%(directories[k],prefixes[k],gradientFiles[j])
            commands += ['./paste_control_forcing %s %s %d %d %d'                                                           \
                         % (globalGradFiles[j], sliceGradFile, totalTimestep, kOffset, Nts)]
    commandString += scriptor.parallelLoopCommand(commands,'paste','paste_control_forcing')

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

def forwardRunCommand(baseDirectory='x0',zeroControlForcing=False):
    bdir = baseDirectory

    commandString = ''
    commandString += distributeCommand(bdir,zeroControlForcing)

    commands, commandDirs = [], []
    for k in range(Nsplit):
        commandDirs += ['%s/%s'%(bdir,directories[k])]
        commands += ['./forward --input %s' % inputFiles[k]]
    commandString += scriptor.parallelLoopCommand(commands,'forward',
                                                  'forward',directories=commandDirs)

    commands = []
    for k in range(Nsplit):
        matchingFile = '%s/%s/%s'%(bdir,directories[k-1],matchingForwardFiles[k-1])
        icFile = '%s/%s'%(bdir,icFiles[k])
        commands += ['./qfile_zaxpy %s %.16E %s %s --input %s'                              \
                    % (diffFiles[k],-1.0,matchingFile,icFile,globalInputFile)]
    commandString += scriptor.parallelLoopCommand(commands,'qfile-zaxpy',
                                                  'forward_qfile_zaxpy')

    commands = []
    for k in range(Nsplit):
        commands += ['./spatial_inner_product %s %s --output %s --input %s'                 \
                    % (diffFiles[k], diffFiles[k], diffOutputFiles[k], globalInputFile)]
    commandString += scriptor.parallelLoopCommand(commands,'qfile-zxdoty',
                                                  'forward_qfile_zxdoty')

    if (useLagrangian):
        commands = []
        for k in range(Nsplit):
            commands += ['./spatial_inner_product %s %s --output %s --input %s'                 \
                        % (diffFiles[k], lagrangianFiles[k], lagrangianOutputFiles[k], globalInputFile)]
        commandString += scriptor.parallelLoopCommand(commands,'qfile-zxdoty',
                                                      'forward_qfile_lagrangian')

    return commandString

def adjointRunCommand(baseDirectory='x0'):
    bdir = baseDirectory

    commandString = ''
    commandString += generalSetOptionCommand

    targetInputFiles = ['%s/%s/%s'%(bdir,dir,file) for dir, file in zip(directories,inputFiles)]

    if (ignoreObjective):
        commands = []
        for k in range(Nsplit):
            commands += ['./setOption.sh %s "adjoint_forcing_switch" "false"' % targetInputFiles[k]]
        commandString += scriptor.nonMPILoopCommand(commands,'magudi_option_disable_adjoint_forcing')

    commands = []
    for k in range(Nsplit):
        matchingAdjointFile = '%s/%s/%s'%(bdir,directories[k-1],matchingAdjointFiles[k-1])
        commands += ['./qfile_zaxpy %s %.16E %s --zero --input %s'                                      \
        % (matchingAdjointFile, -matchingConditionWeight[k], diffFiles[k], globalInputFile)]
    commandString += scriptor.parallelLoopCommand(commands,'qfile-zaxpy',
                                            'adjoint_run_qfile_zaxpy')

    if (useLagrangian):
        commands = []
        for k in range(Nsplit):
            matchingAdjointFile = '%s/%s/%s'%(bdir,directories[k-1],matchingAdjointFiles[k-1])
            commands += ['./qfile_zaxpy %s %.16E %s %s --input %s'                                      \
            % (matchingAdjointFile, 1.0, lagrangianFiles[k], matchingAdjointFile, globalInputFile)]
        commandString += scriptor.parallelLoopCommand(commands,'qfile-zaxpy',
                                                'adjoint_run_qfile_lagrangian')

    commands = []
    for k in range(Nsplit):
        commands += ['./setOption.sh %s "adjoint_nonzero_initial_condition" "true"'                 \
                    % targetInputFiles[k]]
    commandString += scriptor.nonMPILoopCommand(commands,'magudi_option_nonzero_initial_condition_true')

    commands, commandDirs = [], []
    for k in range(Nsplit):
        commandDirs += ['%s/%s' % (bdir,directories[k])]
        commands += ['./adjoint --input %s' % inputFiles[k]]
    commandString += scriptor.parallelLoopCommand(commands,'adjoint',
                                                'adjoint',directories=commandDirs)

    # TODO: Either check the last adjoint run, or complete the periodic optimization!!
    if (ignoreObjective):
        commands = []
        for k in range(Nsplit):
            commands += ['./setOption.sh %s "adjoint_forcing_switch" "true"' % targetInputFiles[k]]
        commandString += scriptor.nonMPILoopCommand(commands,'magudi_option_enable_adjoint_forcing')

    commands = []
    for k in range(Nsplit):
        commands += ['./setOption.sh %s "adjoint_nonzero_initial_condition" "false"'                            \
                    % targetInputFiles[k]]
    commandString += scriptor.nonMPILoopCommand(commands,'magudi_option_nonzero_initial_condition_false')

    commands = []
    for k in range(Nsplit):
        idx = NcontrolRegion + k
        icAdjointFile = '%s/%s/%s' % (bdir,directories[k],icAdjointFiles[k])
        commands += ['./qfile_zaxpy %s %.16E %s %s --input %s'                                      \
                     % (globalGradFiles[idx], matchingConditionWeight[k], diffFiles[k], icAdjointFile, globalInputFile)]
    commandString += scriptor.parallelLoopCommand(commands,'qfile-zaxpy',
                                            'adjoint_run_ic')

    if (useLagrangian):
        commands = []
        for k in range(Nsplit):
            idx = NcontrolRegion + k
            commands += ['./qfile_zaxpy %s %.16E %s %s --input %s'                                     \
                         % (globalGradFiles[idx], -1.0, lagrangianFiles[k], globalGradFiles[idx], globalInputFile)]
        commandString += scriptor.parallelLoopCommand(commands,'qfile-zaxpy',
                                                'adjoint_run_ic_lagrangian')

    return commandString

def dggCommand():
    commandString = ''

    commands = []
    for j in range(NcontrolRegion):
        commands += ['./zwxmwy %s %s %s %s %s'                                                    \
        %(dggFiles[j], globalGradFiles[j], globalGradFiles[j], previousGradFiles[j], globalNormFiles[j])]
    commandString += scriptor.parallelLoopCommand(commands,'zxdoty',
                                            'dgg_zwxmwy')

    commands = []
    for k in range(Nsplit):
        idx = NcontrolRegion + k
        diffFile = '%s/%s.diff.adjoint.q'%(DIFFDIR,prefixes[k])
        prevICGradFile = '%s' % (previousGradFiles[idx])
        icGradFile = '%s' % (globalGradFiles[idx])
        commands += ['./qfile_zaxpy %s %.16E %s %s --input %s'                                    \
                            % (diffFile, -1.0, prevICGradFile, icGradFile, globalInputFile)]
    commandString += scriptor.parallelLoopCommand(commands,'qfile-zaxpy',
                                            'dgg_qfile_zaxpy')

    commands = []
    for k in range(Nsplit):
        idx = NcontrolRegion + k
        diffFile = '%s/%s.diff.adjoint.q'%(DIFFDIR,prefixes[k])
        commands += ['./spatial_inner_product %s %s --output %s --input %s'                                \
                            % (diffFile, globalGradFiles[idx], dggFiles[idx], globalInputFile)]
    commandString += scriptor.parallelLoopCommand(commands,'qfile-zxdoty',
                                            'dgg_qfile_zxdoty')

    return commandString

def updateLagrangian(weight,initial=False):
    commandString = ''

    commands = []
    for k in range(Nsplit):
        a = weight[k]
        command = './qfile_zaxpy ' + lagrangianFiles[k] + ' ' + "{:.16E}".format(-a) + ' ' + diffFiles[k]
        if (initial):
            command += ' --zero'
        else:
            command += ' '+lagrangianFiles[k]
        commands += [command]
    commandString += bashParallelLoopCommand(commands,'qfile-zaxpy',
                                            'update_lagrangian')

    return commandString
