from base import *

def beforeLinmin(initial, zeroBaseline):
    import pandas as pd
    import subprocess

    J0, subJ0 = QoI()
    gg, subgg = readInnerProduct(ggFiles)

    data = [[J0]+list(subJ0)]
    J_new_df = pd.DataFrame(data,columns=forwardLogColumns)
    data = [[gg]+list(subgg)]
    dJ_new_df = pd.DataFrame(data,columns=gradientLogColumns)

    if (initial):
        J_new_df.to_csv(forwardLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)
        dJ_new_df.to_csv(gradientLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)

        df = pd.DataFrame([np.nan],columns=['reduction'])
        df.to_csv(CGLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)

        commands = []
        for k in range(NcontrolSpace):
            commands += ['cp %s %s' % (globalGradFiles[k], globalConjugateGradientFiles[k])]
        commandString = bashParallelCopyCommand(commands,'initial-conjugate-gradient')
        fID = open(globalCommandFile,'w')
        fID.write(commandString)
        fID.close()

        command = 'python3 '+decisionMaker+' 2'
        if(zeroBaseline):
            command += ' -zero_baseline'
        fID = open(decisionMakerCommandFile,'w')
        fID.write(command+'\n')
        fID.close()

        print ('Initial line minimization is ready. Run mnbrak and linmin procedures.')
        return 0

    df = pd.read_csv(forwardLog, sep='\t', header=0)
    J1 = df.at[df.index[-1],'total']
    df = df.append(J_new_df)
    df.to_csv(forwardLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)

    df = pd.read_csv(gradientLog, sep='\t', header=0)
    gg1 = df.at[df.index[-1],'total']
    gg0 = df.at[df.index[0],'total']
    df = df.append(dJ_new_df)
    df.to_csv(gradientLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)

    df = pd.read_csv(CGLog, sep='\t', header=0)
    df_addendum = pd.DataFrame([(1.0-J0/J1)], columns=['reduction'])
    df = df.append(df_addendum)
    df.to_csv(CGLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)

    #reduction = gg/gg0
    #if (reduction<=tol):
    if (gg<=tol):
        print ('FRPRMN - after linmin: conjugate-gradient optimization is finished.')
        commandFile = open(globalCommandFile,'w')
        commandFile.write('exit 1\n')
        commandFile.close()
        commandFile = open(decisionMakerCommandFile,'w')
        commandFile.write('exit 1\n')
        commandFile.close()
        return 0

    dgg, dummy = readInnerProduct(dggFiles)

    # Fletcher-Reeves
    gamma1 = gg/gg1
    # Polak-Ribiere
    gamma = dgg/gg1
    if (gamma > gamma1):
        gamma = gamma1
    elif (gamma < -gamma1):
        gamma = -gamma1

    commandString = zaxpyCommand(globalConjugateGradientFiles,gamma,previousCGFiles,globalGradFiles)
    commandString += '\n'
    commandFile = open(globalCommandFile,'w')
    commandFile.write(commandString)
    commandFile.close()

    commandFile = open(decisionMakerCommandFile,'w')
    command = 'python3 '+decisionMaker+' 2'
    if(zeroBaseline):
        command += ' -zero_baseline'
    commandFile.write(command+'\n')
    commandFile.close()

    print ('line minimization is ready. Run mnbrak and linmin procedures.')
    return 0

def afterLinmin(zeroBaseline):
    import pandas as pd
    import subprocess

    #copy line minimization log
    import os
    numFiles = len(os.listdir('./linminLog/'))
    subprocess.check_call('mkdir -p linminLog/%d'%numFiles,shell=True)
    subprocess.check_call('cp %s linminLog/%d/'%(lineMinLog,numFiles), shell=True)
    for k in range(Nsplit):
        costFunctionalFile = 'b/%s/%s.cost_functional.txt'%(directories[k],prefixes[k])
        subprocess.check_call('cp %s linminLog/%d/'%(costFunctionalFile,numFiles), shell=True)
        if (k<Nsplit-1):
            subprocess.check_call('cp %s linminLog/%d/'%(diffOutputFiles[k],numFiles), shell=True)

    commandString = ''
    for k in range(NcontrolSpace):
        commandString += 'cp b/%s x0/ \n' % (globalControlSpaceFiles[k])
        commandString += 'mv %s %s \n' % (globalGradFiles[k],previousGradFiles[k])
        commandString += 'mv %s %s \n' % (globalConjugateGradientFiles[k],previousCGFiles[k])
    commandString += '\n'
    if( zeroBaseline ):
        commandString += generalSetOptionCommand
        targetInputFiles = ['x0/%s/%s'%(dir,file) for dir, file in zip(directories,inputFiles)]
        commands = []
        for k in range(Nsplit):
            commands += ['setOption %s "controller_switch" true' % targetInputFiles[k]]
        commandString += bashParallelCopyCommand(commands,'magudi_option_turn_on_controller')

    commandString += forwardRunCommand()
    commandString += '\n'
    commandString += adjointRunCommand()
    commandString += '\n'
    commandString += gatherControlForcingGradientCommand()
    commandString += '\n'
    commandString += innerProductCommand(globalGradFiles,globalGradFiles,ggFiles)

    commands = []
    for k in range(Nsplit):
        costSensitivityFile = 'x0/%s/%s.cost_sensitivity.txt'%(directories[k],prefixes[k])
        commands += ['cp %s linminLog/%d/ '%(costSensitivityFile,numFiles)]
    commandString += bashParallelCopyCommand(commands,'copy_gradient_history')

    # Polak-Ribiere
    commandString += dggCommand()
    fID = open(globalCommandFile,'w')
    fID.write(commandString)
    fID.close()
    commandFile = open(decisionMakerCommandFile,'w')
    commandFile.write('python3 '+decisionMaker+' 1 \n')
    commandFile.close()

    print ('FRPRMN - after linmin: postprocessing is finished. Run new forward/adjoint simulations.')
    return 1
