from base import *
from filenames_newton import *

import numpy as np
import subprocess
import pandas as pd

def forward(zeroControl=False):
    bdir = 'x0'

    commandString = ''
    commandString += forwardRunCommand(bdir,zeroControl)
    fID = open(globalCommandFile,'w')
    fID.write(commandString)
    fID.close()

    commandFile = open(decisionMakerCommandFile,'w')
    command = 'python3 %s residual' % decisionMaker
    if(zeroControl):
        command += ' -zero_control'
    commandFile.write(command+'\n')
    commandFile.close()

    print ('forward procedure is ready.')
    return True

def residual(zeroControl=False):

    dummy, subrr = QoI()

    # Update or create optimization log.
    data = [[dummy]+list(subrr)]
    J_new_df = pd.DataFrame(data,columns=forwardLogColumns)
    from os import path
    if (path.exists(forwardLog)):
        df = pd.read_csv(forwardLog, sep='\t', header=0)
        J_new_df = df.append(J_new_df)
    J_new_df.to_csv(forwardLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)

    # Check the residual.
    rr = np.zeros(Nsplit+1)
    rr[1:] = np.roll( 2.0 * subrr[Nsplit:2*Nsplit], 1 )
    if (not periodicSolution):
        rr[1] = 0.0
    rr[0] = np.sum( rr[1:] )

    if (rr[0] < eps):
        print ('Residual is below tolerance: %.5E' % rr[0])
        commandFile = open(globalCommandFile,'w')
        commandFile.write('exit 1\n')
        commandFile.close()
        commandFile = open(decisionMakerCommandFile,'w')
        commandFile.write('exit 1\n')
        commandFile.close()
        return True

    # Create a new residual log.
    res_new_df = pd.DataFrame([list(rr)],columns=residualLogColumns)
    res_new_df.to_csv(residualLog, float_format='%.16E', encoding='utf-8',
                                                sep='\t', mode='w', index=False)

    commandString = ''

    # Copy forward run files
    import os
    numFiles = len(os.listdir('./linminLog/'))
    subprocess.check_call('mkdir -p linminLog/%d'%numFiles,shell=True)
    commands = []
    for k in range(Nsplit):
        costFunctionalFile = 'x0/%s/%s.cost_functional.txt'%(directories[k],prefixes[k])
        commands += ['cp %s linminLog/%d/'%(costFunctionalFile,numFiles)]
        commands += ['cp %s linminLog/%d/'%(diffOutputFiles[k],numFiles)]
    commandString += bashParallelCopyCommand(commands,'saving_residual_files')

    # move the residual files.
    # diffFiles[k] = Q_{k+1}^{-} - Q_{k+1}^{+}
    # diff file is negative of the residual in the formulation.
    commands = []
    for k in range(Nsplit):
        if ((k==0) and (not periodicSolution)):
            continue
        commands += ['cp %s %s' % (diffFiles[k-1], residualFiles[k])]
    commandString += bashParallelCopyCommand(commands,'copy_residual')

    if (not periodicSolution):
        commands = []
        commands += ['./qfile_zaxpy %s %.16E %s --zero' % (residualFiles[0],0.0,diffFiles[-1])]
        commandString += bashSerialLoopCommand(commands,'qfile-zaxpy',
                                                prefix='zero_residual')

    commands = []
    for k in range(Nsplit):
        idx = NcontrolRegion + k
        commands += ['cp %s %s'                                                 \
                     % (residualFiles[k],globalConjugateGradientFiles[idx])]
    commandString += bashParallelCopyCommand(commands,'copy_conjugate_gradient')

    commands = []
    for k in range(Nsplit):
        commands += ['./qfile_zaxpy %s %.16E %s --zero'                                \
                    % (newtonAdjointFiles[k],0.0,residualFiles[k])]
    commandString += bashParallelLoopCommand(commands,'qfile-zaxpy',
                                            'zero_newton_adjoint')

    fID = open(globalCommandFile,'w')
    fID.write(commandString)
    fID.close()

    commandFile = open(decisionMakerCommandFile,'w')
    command = 'python3 %s adjoint' % decisionMaker
    if(zeroControl):
        command += ' -zero_control'
    commandFile.write(command+'\n')
    commandFile.close()

    print ('forward procedure is finished. residual procedure is ready.')
    return True

def adjoint(zeroControl=False):
    bdir = 'x0'

    commandString = ''

    commands = []
    for k in range(Nsplit):
        if ((k==0) and (not periodicSolution)):
            continue
        idx = NcontrolRegion + k
        matchingAdjointFile = '%s/%s/%s'                                        \
        % (bdir,directories[k-1],matchingAdjointFiles[k-1])
        commands += ['cp %s %s'                                                 \
                    % (globalConjugateGradientFiles[idx],matchingAdjointFile)]
    commandString += bashParallelCopyCommand(commands,'adjoint_final_condition')

    commandString += generalSetOptionCommand
    targetInputFiles = ['%s/%s/%s'%(bdir,dir,file) for dir, file in zip(directories,inputFiles)]

    commands = []
    for k in range(Nsplit):
        if ((k==Nsplit-1) and (not periodicSolution)):
            continue
        commands += ['setOption %s "adjoint_nonzero_initial_condition" "true"'                 \
                    % targetInputFiles[k]]
    commandString += bashParallelCopyCommand(commands,'magudi_option_nonzero_initial_condition_true')

    commands = []
    for k in range(Nsplit):
        if ((k==Nsplit-1) and (not periodicSolution)):
            continue
        commands += ['setOption %s "adjoint_forcing_switch" "false"'                                    \
                    % targetInputFiles[k]]
    commandString += bashParallelCopyCommand(commands,'magudi_option_disable_adjoint_forcing')

    commandString += '\n'

    commands, commandDirs = [], []
    for k in range(Nsplit):
        if ((k==Nsplit-1) and (not periodicSolution)):
            continue
        commandDirs += ['%s/%s' % (bdir,directories[k])]
        commands += ['./adjoint --input %s' % inputFiles[k]]
    commandString += bashParallelLoopCommand(commands,'adjoint',
                                            'adjoint_run',directories=commandDirs)

    commands = []
    for k in range(Nsplit):
        if ((k==Nsplit-1) and (not periodicSolution)):
            continue
        commands += ['setOption %s "adjoint_nonzero_initial_condition" "false"'                 \
                    % targetInputFiles[k]]
    commandString += bashParallelCopyCommand(commands,'magudi_option_nonzero_initial_condition_false')

    commands = []
    for k in range(Nsplit):
        if ((k==Nsplit-1) and (not periodicSolution)):
            continue
        commands += ['setOption %s "adjoint_forcing_switch" "true"'                                    \
                    % targetInputFiles[k]]
    commandString += bashParallelCopyCommand(commands,'magudi_option_enable_adjoint_forcing')

    commandString += '\n'

    # NOTE: output variables have the opposite sign of the formulation.
    commands = []
    for k in range(Nsplit):
        idx = NcontrolRegion + k
        icLinearizedFile = '%s/%s/%s' % (bdir,directories[k],icLinearizedFiles[k])
        if ((k==Nsplit-1) and (not periodicSolution)):
            icAdjointFile = '--zero'
        else:
            icAdjointFile = '%s/%s/%s' % (bdir,directories[k],icAdjointFiles[k])
        commands += ['./qfile_zaxpy %s %.16E %s %s --input %s'                                      \
        % (icLinearizedFile, -1.0, globalConjugateGradientFiles[idx], icAdjointFile, globalInputFile)]
    commandString += bashParallelLoopCommand(commands,'qfile-zaxpy','linearized_run_ic')

    fID = open(globalCommandFile,'w')
    fID.write(commandString)
    fID.close()

    commandFile = open(decisionMakerCommandFile,'w')
    command = 'python3 %s linearized' % decisionMaker
    if(zeroControl):
        command += ' -zero_control'
    commandFile.write(command+'\n')
    commandFile.close()

    print ('residual procedure is finished. adjoint procedure is ready.')
    return True

def linearized(zeroControl=False):
    bdir = 'x0'

    commandString = ''

    if (not periodicSolution):
        commands = []
        icLinearizedFile = '%s/%s/%s' % (bdir,directories[0],icLinearizedFiles[0])
        commands += ['./qfile_zaxpy %s %.16E %s --zero --input %s'                              \
        % (icLinearizedFile, 0.0, icLinearizedFile, globalInputFile)]
        commandString += bashSerialLoopCommand(commands,'qfile-zaxpy',prefix='zero_linearized_ic')

    commands, commandDirs = [], []
    for k in range(Nsplit):
        if ((k==Nsplit-1) and (not periodicSolution)):
            continue
        commandDirs += ['%s/%s' % (bdir,directories[k])]
        commands += ['./linearized --input %s' % inputFiles[k]]
    commandString += bashParallelLoopCommand(commands,'adjoint','linearized',directories=commandDirs)

    # NOTE: output linearized variables have the opposite sign of the formulation.
    commands = []
    for k in range(Nsplit):
        icLinearizedFile = '%s/%s/%s' % (bdir,directories[k],icLinearizedFiles[k])
        if ((k==0) and (not periodicSolution)):
            amp = 0.0
            matchingLinearizedFile = '--zero'
        else:
            amp = - 1.0
            matchingLinearizedFile = '%s/%s/%s'                                     \
            % (bdir,directories[k-1],matchingLinearizedFiles[k-1])
        commands += ['./qfile_zaxpy %s %.16E %s %s --input %s'                                      \
        % (diffResidualFiles[k], amp, icLinearizedFile, matchingLinearizedFile, globalInputFile)]
    commandString += bashParallelLoopCommand(commands,'qfile-zaxpy','linearized_run_diff')

    # inner product for line minimization step
    commands = []
    for k in range(Nsplit):
        idx = NcontrolRegion + k
        commands += ['./spatial_inner_product %s %s --output %s --input %s'                     \
        % (globalConjugateGradientFiles[idx],diffResidualFiles[k],dggFiles[idx],globalInputFile)]
    commandString += bashParallelLoopCommand(commands,'qfile-zxdoty','linearized_run_inner_product')

    fID = open(globalCommandFile,'w')
    fID.write(commandString)
    fID.close()

    commandFile = open(decisionMakerCommandFile,'w')
    command = 'python3 %s linmin' % decisionMaker
    if(zeroControl):
        command += ' -zero_control'
    commandFile.write(command+'\n')
    commandFile.close()

    print ('adjoint procedure is finished. linearized run procedure is ready.')
    return True

def linmin(zeroControl=False):

    dummy, subgg = readInnerProduct(dggFiles)
    pAATp = np.sum( subgg[-Nsplit:] )

    df = pd.read_csv(residualLog, sep='\t', header=0)
    rr = df.at[df.index[-1],'total']
    a = rr / pAATp

    commandString = ''

    commands = []
    for k in range(Nsplit):
        if ((k==0) and (not periodicSolution)):
            continue
        idx = NcontrolRegion + k
        commands += ['./qfile_zaxpy %s %.16E %s %s'                                \
        % (newtonAdjointFiles[k], a ,globalConjugateGradientFiles[idx],newtonAdjointFiles[k])]
    commandString += bashParallelLoopCommand(commands,'qfile-zaxpy','newton_adjoint_update')

    commands = []
    for k in range(Nsplit):
        if ((k==0) and (not periodicSolution)):
            continue
        idx = NcontrolRegion + k
        commands += ['./qfile_zaxpy %s %.16E %s %s'                                \
        % (residualFiles[k], - a ,diffResidualFiles[k],residualFiles[k])]
    commandString += bashParallelLoopCommand(commands,'qfile-zaxpy','residual_update')

    commands = []
    for k in range(Nsplit):
        idx = NcontrolRegion + k
        commands += ['./spatial_inner_product %s %s --output %s --input %s'                     \
        % (residualFiles[k],residualFiles[k],ggFiles[idx],globalInputFile)]
    commandString += bashParallelLoopCommand(commands,'qfile-zxdoty','residual_inner_product')

    fID = open(globalCommandFile,'w')
    fID.write(commandString)
    fID.close()

    commandFile = open(decisionMakerCommandFile,'w')
    command = 'python3 %s cgstep' % decisionMaker
    if(zeroControl):
        command += ' -zero_control'
    commandFile.write(command+'\n')
    commandFile.close()

    print ('linearized run procedure is finished. linmin procedure is ready.')
    return True

def cgstep(zeroControl=False):

    dummy, subgg = readInnerProduct(ggFiles)
    rr = np.zeros(Nsplit+1)
    rr[-Nsplit:] = subgg[-Nsplit:]
    rr[0] = np.sum( rr[-Nsplit:] )

    df = pd.read_csv(residualLog, sep='\t', header=0)
    rr1 = df.at[df.index[-1],'total']
    rr0 = df.at[df.index[0],'total']

    res_new_df = pd.DataFrame([list(rr)],columns=residualLogColumns)
    df = df.append(res_new_df)
    df.to_csv(residualLog, float_format='%.16E', encoding='utf-8',
                                    sep='\t', mode='w', index=False)

    if (rr[0]/rr0 < linminTol):
        print ('inexact newton step is computed. ready to update control parameters.')
        commandFile = open(globalCommandFile,'w')
        commandFile.close()
        commandFile = open(decisionMakerCommandFile,'w')
        command = 'python3 %s newtonstep' % decisionMaker
        if(zeroControl):
            command += ' -zero_control'
        commandFile.write(command+'\n')
        commandFile.close()
        return True

    be = rr[0] / rr1

    commandString = ''

    commands = []
    for k in range(Nsplit):
        if ((k==0) and (not periodicSolution)):
            continue
        idx = NcontrolRegion + k
        commands += ['./qfile_zaxpy %s %.16E %s %s'                                \
        % (globalConjugateGradientFiles[idx],be,globalConjugateGradientFiles[idx],residualFiles[k])]
    commandString += bashParallelLoopCommand(commands,'qfile-zaxpy','conjugate_gradient_update')

    fID = open(globalCommandFile,'w')
    fID.write(commandString)
    fID.close()

    commandFile = open(decisionMakerCommandFile,'w')
    command = 'python3 %s adjoint' % decisionMaker
    if(zeroControl):
        command += ' -zero_control'
    commandFile.write(command+'\n')
    commandFile.close()

    print ('linmin procedure is finished. cgstep procedure is ready.')
    return True

def newtonstep(zeroControl=False):
    bdir = 'x0'

    #copy residual log
    import os
    numFiles = len(os.listdir('./linminLog/'))
    subprocess.check_call('cp %s linminLog/%d/'%(residualLog,numFiles-1), shell=True)

    commandString = ''

    commands = []
    for k in range(Nsplit):
        if ((k==0) and (not periodicSolution)):
            continue
        matchingAdjointFile = '%s/%s/%s'                                        \
        % (bdir,directories[k-1],matchingAdjointFiles[k-1])
        commands += ['cp %s %s'                                                 \
                    % (newtonAdjointFiles[k],matchingAdjointFile)]
    commandString += bashParallelCopyCommand(commands,'adjoint_final_condition')

    commandString += generalSetOptionCommand
    targetInputFiles = ['%s/%s/%s'%(bdir,dir,file) for dir, file in zip(directories,inputFiles)]

    commands = []
    for k in range(Nsplit):
        if ((k==Nsplit-1) and (not periodicSolution)):
            continue
        commands += ['setOption %s "adjoint_nonzero_initial_condition" "true"'                 \
                    % targetInputFiles[k]]
    commandString += bashParallelCopyCommand(commands,'magudi_option_nonzero_initial_condition_true')

    commands = []
    for k in range(Nsplit):
        if ((k==Nsplit-1) and (not periodicSolution)):
            continue
        commands += ['setOption %s "adjoint_forcing_switch" "false"'                                    \
                    % targetInputFiles[k]]
    commandString += bashParallelCopyCommand(commands,'magudi_option_disable_adjoint_forcing')

    commandString += '\n'

    commands, commandDirs = [], []
    for k in range(Nsplit):
        if ((k==Nsplit-1) and (not periodicSolution)):
            continue
        commandDirs += ['%s/%s' % (bdir,directories[k])]
        commands += ['./adjoint --input %s' % inputFiles[k]]
    commandString += bashParallelLoopCommand(commands,'adjoint',
                                            'newton-adjoint',directories=commandDirs)

    commands = []
    for k in range(Nsplit):
        if ((k==Nsplit-1) and (not periodicSolution)):
            continue
        commands += ['setOption %s "adjoint_nonzero_initial_condition" "false"'                 \
                    % targetInputFiles[k]]
    commandString += bashParallelCopyCommand(commands,'magudi_option_nonzero_initial_condition_false')

    commands = []
    for k in range(Nsplit):
        if ((k==Nsplit-1) and (not periodicSolution)):
            continue
        commands += ['setOption %s "adjoint_forcing_switch" "true"'                                    \
                    % targetInputFiles[k]]
    commandString += bashParallelCopyCommand(commands,'magudi_option_enable_adjoint_forcing')

    commandString += '\n'

    # NOTE: output variables have the opposite sign of the formulation.
    commands = []
    for k in range(Nsplit):
        idx = NcontrolRegion + k
        if ((k==Nsplit-1) and (not periodicSolution)):
            icAdjointFile = '--zero'
        else:
            icAdjointFile = '%s/%s/%s' % (bdir,directories[k],icAdjointFiles[k])
        commands += ['./qfile_zaxpy %s %.16E %s %s --input %s'                                      \
        % (newtonFiles[idx], -1.0, newtonAdjointFiles[k], icAdjointFile, globalInputFile)]
    commandString += bashParallelLoopCommand(commands,'qfile-zaxpy','newton_step_ic')

    if (not periodicSolution):
        commands = []
        idx = NcontrolRegion
        commands += ['./qfile_zaxpy %s %.16E %s --zero --input %s'                              \
        % (newtonFiles[idx], 0.0, newtonFiles[idx], globalInputFile)]
        commandString += bashSerialLoopCommand(commands,'qfile-zaxpy',prefix='zero_newton_ic')
        if (not zeroControl):
            commands = []
            for j in range(NcontrolRegion):
                zeroSlice = 'x0/%s/%s%s'                                        \
                % (directories[-1],prefixes[-1],gradientFiles[j])
                tempSlice = 'x0/%s/%s%s'                                        \
                % (directories[0],prefixes[0],gradientFiles[j])
                commands += ['./zaxpy %s 0.0 %s' % (zeroSlice,tempSlice)]
            commandString += bashParallelLoopCommand(commands,'zaxpy',
                                        'newton_step_control_forcing')

    if (not zeroControl):
        commandString += gatherControlForcingGradientCommand()
        commands = []
        for j in range(NcontrolRegion):
            commands += ['mv %s %s' % (globalGradFiles[j],newtonFiles[j])]
        commandString += bashParallelCopyCommand(commands,'rename_newton_step_control_forcing')

    commands = []
    target = ['x0/'+file for file in globalControlSpaceFiles]
    for k in range(Nsplit):
        if ((k==0) and (not periodicSolution)):
            continue
        idx = NcontrolRegion + k
        commands += ['./qfile_zaxpy %s %.16E %s %s --input %s'                                      \
        % (target[idx], -1.0, newtonFiles[idx], target[idx], globalInputFile)]
    commandString += bashParallelLoopCommand(commands,'qfile-zaxpy','newton_update_ic')

    if (not zeroControl):
        commands = []
        for j in range(NcontrolRegion):
            commands += ['./zaxpy %s -1.0 %s %s' % (target[j],newtonFiles[j],target[j])]
        commandString += bashParallelLoopCommand(commands,'zaxpy',
                                                'newton_update_control_forcing')

    fID = open(globalCommandFile,'w')
    fID.write(commandString)
    fID.close()

    commandFile = open(decisionMakerCommandFile,'w')
    command = 'python3 %s forward' % decisionMaker
    if(zeroControl):
        command += ' -zero_control'
    commandFile.write(command+'\n')
    commandFile.close()
    return True
