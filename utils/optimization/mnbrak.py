from base import *

def setupInitialSteps(forwardFilename, CGFilenames, controlForcingFilenames, zeroBaseline=True):
    from base import readScalar
    import subprocess
    import pandas as pd
    NumCGFile = len(CGFilenames)
    if (len(controlForcingFilenames)!=NumCGFile):
        print ('MNBRAK - setup initial steps: numbers of control space files do not match!')
        return -1

    J0 = readScalar(forwardFilename)
    subprocess.check_call('cp '+forwardFilename+' '+str(NumSearch+1)+'/',shell=True)
    if (not zeroBaseline):
        for k in range(NumCGFile):
            subprocess.check_call('cp '+controlForcingFilenames[k]+' '+str(NumSearch+1)+'/ &',shell=True)

    steps, Js = np.zeros(NumSearch+1), np.zeros(NumSearch+1)
    steps[1:NumSearch+1] = initial_step*golden_ratio**np.linspace(0.0,NumSearch-1,NumSearch)
    Js[0] = J0

    commandFile = open(commandFilename,'w')
    for k in range(NumSearch):
        for i in range(NumCGFile):
            command = 'srun -n '+str(NumProcs)+' ./zaxpy '                                                 \
                        +str(k+1)+'/'+controlForcingFilenames[i]+" "                                    \
                        +"{:.16E}".format(-steps[k+1])+" "+CGFilenames[i]
            if (not zeroBaseline):
                command += " "+controlForcingFilenames[i]
            commandFile.write(command+'\n')
    commandFile.write('bash intermediate_forward_runs.sh\n')
    commandFile.close()
    commandFile = open(decisionMakerCommandFilename,'w')
    command = 'python '+decisionMaker+' 3'
    if(zeroBaseline):
        command += ' -zero_baseline'
    commandFile.write(command+'\n')
    commandFile.close()

    data = {'step':steps,'QoI':Js,'directory index':[NumSearch+1]+list(range(1,NumSearch+1))}
    df = pd.DataFrame(data)
    df.to_csv(lineMinLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)
    print ('Initial steps written in command. Run '+commandFilename+'.')
    return 0

def NextMnbrak(forwardFilename, CGFilenames, controlForcingFilenames, zeroBaseline=True, retry=False):
    from base import readScalar
    import subprocess
    import pandas as pd
    # from IPython.display import display, HTML
    NumCGFile = len(CGFilenames)
    if (len(controlForcingFilenames)!=NumCGFile):
        print ('MNBRAK - next mnbrak: numbers of control space files do not match!')
        return -1

    if (retry):
        df = pd.read_csv(lineMinLog, sep='\t', header=0)
        commandFile = open(commandFilename,'w')
        for k in range(NumSearch):
            stepK = float(df.loc[df['directory index']==k+1,'step'])
            for i in range(NumCGFile):
                command = 'srun -n '+str(NumProcs)+' ./zaxpy '                                                 \
                            +str(k+1)+'/'+controlForcingFilenames[i]+' '                                    \
                            +"{:.16E}".format(-stepK)+' '+CGFilenames[i]
                if (not zeroBaseline):
                    command += ' '+controlForcingFilenames[i]
                commandFile.write(command+'\n')
        commandFile.write('bash intermediate_forward_runs.sh\n')
        commandFile.close()
        commandFile = open(decisionMakerCommandFilename,'w')
        command = 'python '+decisionMaker+' 3'
        if(zeroBaseline):
            command += ' -zero_baseline'
        commandFile.write(command+'\n')
        commandFile.close()

        print ('submitted all zaxpy works. Wait until they finish.')
        print (df[df['directory index']>NumSearch])
        print ('MNBRAK: rerunning zaxpy works. Run '+commandFilename+'.')
        return 3
    

    df = collectQoIs(lineMinLog, forwardFilename, NumSearch)
    steps = np.array(df['step'][df['directory index']>0])
    Js = np.array(df['QoI'][df['directory index']>0])
    dirIdx = np.array(df['directory index'][df['directory index']>0])
    N = steps.size
    minIdx = np.argmin(Js)

    if ( (minIdx==0) | (Js[1]>Js[0]) ):
        bracketSize = (steps[1]-steps[0])/(steps[1]+eps)
        if (bracketSize<tol):
            for k in range(3):
                switchDirectory(dirIdx[k],NumSearch+1+k,df)
                dirIdx = np.array(df['directory index'][df['directory index']>0])
            df.loc[df['directory index']<=NumSearch,'directory index'] = 0
            df.to_csv(lineMinLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)
            print (np.array(df[df['directory index']>NumSearch]))
            print ('MNBRAK: initial mininum bracket is prepared.')

            commandFile = open(commandFilename,'w')
            commandFile.close()
            commandFile = open(decisionMakerCommandFilename,'w')
            command = 'python '+decisionMaker+' 4 -linmin_initial'
            if(zeroBaseline):
                command += ' -zero_baseline'
            commandFile.write(command+'\n')
            commandFile.close()
            return 0
        dirIdx = switchDirectory(1,NumSearch+3,df)
        df.loc[df['directory index']<=NumSearch,'directory index'] = 0

        h = (steps[1]-steps[0])/(NumSearch+1)
        steps = np.linspace(steps[0]+h,steps[1]-h,NumSearch)
        data = {'step':steps, 'QoI':np.ones(NumSearch)*np.nan, 'directory index':list(range(1,NumSearch+1))}
        new_df = pd.DataFrame(data)
        df = df.append(new_df, ignore_index=True)

        commandFile = open(commandFilename,'w')
        for k in range(NumSearch):
            for i in range(NumCGFile):
                command = 'srun -n '+str(NumProcs)+' ./zaxpy '                                                 \
                            +str(k+1)+'/'+controlForcingFilenames[i]+' '                                    \
                            +"{:.16E}".format(-steps[k])+' '+CGFilenames[i]
                if (not zeroBaseline):
                    command += ' '+controlForcingFilenames[i]
                commandFile.write(command+'\n')
        commandFile.write('bash intermediate_forward_runs.sh\n')
        commandFile.close()
        commandFile = open(decisionMakerCommandFilename,'w')
        command = 'python '+decisionMaker+' 3'
        if(zeroBaseline):
            command += ' -zero_baseline'
        commandFile.write(command+'\n')
        commandFile.close()

        df.to_csv(lineMinLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)
        print ('submitted all zaxpy works. Wait until they finish.')
        print (df[df['directory index']>NumSearch])
        print ('MNBRAK: narrowing the bracket - Run intermediate forward simulations.')
        return 1

    elif (minIdx==N-1):
        for k in range(2):
            switchDirectory(dirIdx[minIdx-1+k],NumSearch+1+k,df)
            dirIdx = np.array(df['directory index'][df['directory index']>0])
        df.loc[df['directory index']<=NumSearch,'directory index'] = 0
        df.loc[df['directory index']==NumSearch+3,'directory index'] = 0

        steps = steps[N-1]*golden_ratio**np.linspace(1.0,NumSearch,NumSearch)
        data = {'step':steps, 'QoI':np.ones(NumSearch)*np.nan, 'directory index':list(range(1,NumSearch+1))}
        new_df = pd.DataFrame(data)
        df = df.append(new_df, ignore_index=True)

        commandFile = open(commandFilename,'w')
        for k in range(NumSearch):
            for i in range(NumCGFile):
                command = 'srun -n '+str(NumProcs)+' ./zaxpy '                                                 \
                            +str(k+1)+'/'+controlForcingFilenames[i]+' '                                    \
                            +"{:.16E}".format(-steps[k])+' '+CGFilenames[i]
                if (not zeroBaseline):
                    command += ' '+controlForcingFilenames[i]
                commandFile.write(command+'\n')
        commandFile.write('bash intermediate_forward_runs.sh\n')
        commandFile.close()
        commandFile = open(decisionMakerCommandFilename,'w')
        command = 'python '+decisionMaker+' 3'
        if(zeroBaseline):
            command += ' -zero_baseline'
        commandFile.write(command+'\n')
        commandFile.close()

        df.to_csv(lineMinLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)
        print ('submitted all zaxpy works. Wait until they finish.')
        print (df[df['directory index']>NumSearch])
        print ('MNBRAK: expanding the bracket - Run intermediate forward simulations.')
        return 2

    else:
        for k in range(3):
            switchDirectory(dirIdx[minIdx-1+k],NumSearch+1+k,df)
            dirIdx = np.array(df['directory index'][df['directory index']>0])
#             display(df)
        df.loc[df['directory index']<=NumSearch,'directory index'] = 0
        df.to_csv(lineMinLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)
        print (df[df['directory index']>NumSearch])
        print ('MNBRAK: initial mininum bracket is prepared.')

        commandFile = open(commandFilename,'w')
        commandFile.close()
        commandFile = open(decisionMakerCommandFilename,'w')
        command = 'python '+decisionMaker+' 4 -linmin_initial'
        if(zeroBaseline):
            command += ' -zero_baseline'
        commandFile.write(command+'\n')
        commandFile.close()
        return 0
