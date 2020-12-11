from base import *

def setupInitialSteps(zeroBaseline=True):
    import subprocess
    import pandas as pd

    commandString = ''
    J0, dummy = QoI()

    commands = []
    for k in range(1,Nsplit):
        commands += ['cp x0/%s ./a/ ' % icFiles[k]]
    if (not zeroBaseline):
        for j in range(NcontrolRegion):
            commands += ['cp x0/%s ./a/ ' % globalControlSpaceFiles[j]]
    commandString += scriptor.nonMPILoopCommand(commands,'copy_control_params')

    import os
    numFiles = len(os.listdir('./linminLog/'))
    if (numFiles>0):
        df_temp = pd.read_csv('linminLog/%d/%s'%(numFiles-1,lineMinLog), sep='\t', header=0)
        lastStep = np.array(df_temp['step'][df_temp['directory index']=='b'])[0]

    steps, Js = np.zeros(2), np.zeros(2)
    steps[1] = initial_step if numFiles==0 else lastStep
    Js[0] = J0

    temp = ['x0/'+file for file in globalControlSpaceFiles]
    if (zeroBaseline):
        for j in range(NcontrolRegion):
            temp[j] = ''
    target = globalControlSpaceFiles.copy()
    for k, file in enumerate(target):
        target[k] = 'x/' + file
    commandString += zaxpyCommand(target, -steps[1], globalConjugateGradientFiles, temp)
    commandString += '\n'
    commandString += forwardRunCommand('x')
    commandString += scriptor.purgeDirectoryCommand('x')
    fID = open(globalCommandFile,'w')
    fID.write(commandString)
    fID.close()

    commandFile = open(decisionMakerCommandFile,'w')
    command = 'python3 '+decisionMaker+' 3'
    if(zeroBaseline):
        command += ' -zero_baseline'
    commandFile.write(command+'\n')
    commandFile.close()

    data = {'step':steps,'QoI':Js,'directory index':['a','x']}
    df = pd.DataFrame(data)
    df.to_csv(lineMinLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)
    print ('Initial steps written in command. Run '+globalCommandFile+'.')
    return 0

def NextMnbrak(zeroBaseline=True):
    import subprocess
    import pandas as pd

    df = collectQoI(lineMinLog)
    steps = np.array(df['step'][df['directory index']!='0'])
    Js = np.array(df['step'][df['directory index']!='0'])
    dirIdx = np.array(df['directory index'][df['directory index']!='0'])
    N = steps.size

    a = np.array(df['step'][df['directory index']=='a'])[0]
    Ja = np.array(df['QoI'][df['directory index']=='a'])[0]
    x = np.array(df['step'][df['directory index']=='x'])[0]
    Jx = np.array(df['QoI'][df['directory index']=='x'])[0]
    bExists = sum(dirIdx=='b')
    cExists = sum(dirIdx=='c')

    if (Jx>huge):
        print ('Forward run crashed. Retrying with a smaller step.')
        df.loc[df['directory index']=='x','directory index'] = '0'
        Cr = 1.0 - 1.0/golden_ratio
        if (bExists):
            b = np.array(df['step'][df['directory index']=='b'])[0]
            new_x = b + Cr * (x-b)
        else:
            new_x = a + Cr * (x-a)
    elif (bExists):
        b = np.array(df['step'][df['directory index']=='b'])[0]
        Jb = np.array(df['QoI'][df['directory index']=='b'])[0]

        if (Jx>Jb):
            switchDirectory('x','c',df)

            commandFile = open(globalCommandFile,'w')
            commandFile.write('exit 0\n')
            commandFile.close()
            commandFile = open(decisionMakerCommandFile,'w')
            command = 'python3 '+decisionMaker+' 4 -linmin_initial'
            if(zeroBaseline):
                command += ' -zero_baseline'
            commandFile.write(command+'\n')
            commandFile.close()

            df.to_csv(lineMinLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)
            print (np.array(df[df['directory index']!='0']))
            print ('MNBRAK: initial mininum bracket is prepared.')
            return 0
        else:
            switchDirectory('x','b',df)
            switchDirectory('a','x',df)
            df.loc[df['directory index']=='x','directory index'] = '0'

            if (len(df[(df['QoI']>huge)])>0):
                maxX = np.min(np.array(df['step'][df['QoI']>huge]))
                Cr = 1.0 - 1.0/golden_ratio
                new_x = x + Cr * (maxX-x)
            else:
                new_x = x * golden_ratio if (x<safe_zone) else x * 1.05
            print ('MNBRAK: expanding the bracket - Run intermediate forward simulation.')
    elif (cExists):
        c = np.array(df['step'][df['directory index']=='c'])[0]
        Jc = np.array(df['QoI'][df['directory index']=='c'])[0]

        if (Jx<Ja):
            switchDirectory('x','b',df)

            commandFile = open(globalCommandFile,'w')
            commandFile.write('exit 0\n')
            commandFile.close()
            commandFile = open(decisionMakerCommandFile,'w')
            command = 'python3 '+decisionMaker+' 4 -linmin_initial'
            if(zeroBaseline):
                command += ' -zero_baseline'
            commandFile.write(command+'\n')
            commandFile.close()

            df.to_csv(lineMinLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)
            print (np.array(df[df['directory index']!='0']))
            print ('MNBRAK: initial mininum bracket is prepared.')
            return 0
        else:
            switchDirectory('x','c',df)
            df.loc[df['directory index']=='x','directory index'] = '0'

            new_x = x / golden_ratio
            print ('MNBRAK: narrowing the bracket - Run intermediate forward simulations.')
    else:
        if (Jx>Ja):
            switchDirectory('x','c',df)

            new_x = x / golden_ratio
            print ('MNBRAK: narrowing the bracket - Run intermediate forward simulations.')
        else:
            switchDirectory('x','b',df)

            if (len(df[(df['QoI']>huge)])>0):
                maxX = np.min(np.array(df['step'][df['QoI']>huge]))
                Cr = 1.0 - 1.0/golden_ratio
                new_x = x + Cr * (maxX-x)
            else:
                new_x = x * golden_ratio if (x<safe_zone) else x * 1.05
            print ('MNBRAK: expanding the bracket - Run intermediate forward simulation.')

    data = {'step':[new_x], 'QoI':[np.nan], 'directory index':['x']}
    new_df = pd.DataFrame(data)
    df = df.append(new_df, ignore_index=True)

    commandString = ''
    temp = ['x0/'+file for file in globalControlSpaceFiles]
    if (zeroBaseline):
        for j in range(NcontrolRegion):
            temp[j] = ''
    target = ['x/'+file for file in globalControlSpaceFiles]
    commandString += scriptor.parallelPurgeCommand(target,'purge_target')
    commandString += zaxpyCommand(target, -new_x, globalConjugateGradientFiles, temp)
    commandString += '\n'
    commandString += forwardRunCommand('x')
    commandString += scriptor.purgeDirectoryCommand('x')
    fID = open(globalCommandFile,'w')
    fID.write(commandString)
    fID.close()

    commandFile = open(decisionMakerCommandFile,'w')
    command = 'python3 '+decisionMaker+' 3'
    if(zeroBaseline):
        command += ' -zero_baseline'
    commandFile.write(command+'\n')
    commandFile.close()

    df.to_csv(lineMinLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)
    print (df[df['directory index']!='0'])
    
    return 0

