from base import *

def beforeLinmin(forwardFilename, adjointFilename,
                 gradientFilenames, CGFilenames,
                 normFilenames, initial=True, zeroBaseline=False):
    from base import readScalar
    import pandas as pd
    import subprocess
    NumCGFile = len(CGFilenames)
    if ((len(gradientFilenames)!=NumCGFile)                                                   \
        | (len(normFilenames)!=NumCGFile)):
        print ('FRPRMN - before linmin: numbers of control space files do not match!')
        return -1

    J0, gg = readScalar(forwardFilename), readScalar(adjointFilename)
    if (initial):
        data = [[J0, np.nan, np.nan, gg]]
        df = pd.DataFrame(data,columns=['before linmin','after linmin','reduction','gg'])
        df.to_csv(CGLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)

        commandFile = open(commandFilename,'w')
        for k in range(NumCGFile):
            commandFile.write('cp '+gradientFilenames[k]+' '+CGFilenames[k]+'\n')
        commandFile.close()
        commandFile = open(decisionMakerCommandFilename,'w')
        command = 'python '+decisionMaker+' 2'
        if(zeroBaseline):
            command += ' -zero_baseline'
        commandFile.write(command+'\n')
        commandFile.close()

        print ('Initial line minimization is ready. Run mnbrak and linmin procedures.')
        return 0

    df = pd.read_csv(CGLog, sep='\t', header=0)
    data = {'before linmin':J0, 'after linmin':np.nan, 'reduction':np.nan, 'gg':gg}
    df = df.append(data,ignore_index=True)
    df.to_csv(CGLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)

    # Polak-Ribiere
    dgg = 0.0
    for k in range(NumCGFile):
        command = 'srun -n '+str(NumProcs)+' ./zwxmwy'+' '+dggFilename                                    \
                    +' '+gradientFilenames[k]                                                             \
                    +' '+gradientFilenames[k]+' '+'previous.'+gradientFilenames[k]                        \
                    +' '+normFilenames[k]
        print (command)
        subprocess.check_call(command,shell=True)
        dgg += readScalar(dggFilename)
        print (dgg)
    # # Fletcher-Reeves
    # dgg = df.at[df.index[-1],'gg']
    gg = df.at[df.index[-2],'gg']
    gamma = dgg/gg

    commandFile = open(commandFilename,'w')
    for k in range(NumCGFile):
        command = 'srun -n '+str(NumProcs)+' ./zaxpy '                                              \
                    +CGFilenames[k]+' '                                                             \
                    +"{:.16E}".format(gamma)+' '+'previous.'+CGFilenames[k]+' '+gradientFilenames[k]
        commandFile.write(command+'\n')
    commandFile.close()
    commandFile = open(decisionMakerCommandFilename,'w')
    command = 'python '+decisionMaker+' 2'
    if(zeroBaseline):
        command += ' -zero_baseline'
    commandFile.write(command+'\n')
    commandFile.close()

    print ('line minimization is ready. Run mnbrak and linmin procedures.')
    return 0

def afterLinmin(forwardFilename, adjointFilename,
                 gradientFilenames, CGFilenames,
                 normFilenames, controlForcingFilenames, zeroBaseline=False):
    from base import readScalar
    import pandas as pd
    import subprocess
    NumCGFile = len(CGFilenames)
    if ((len(gradientFilenames)!=NumCGFile)                                                     \
        | (len(normFilenames)!=NumCGFile)                                                       \
        | (len(controlForcingFilenames)!=NumCGFile)):
        print ('FRPRMN - after linmin: numbers of control space files do not match!')
        return -1

    J0 = readScalar(forwardFilename)
    df = pd.read_csv(lineMinLog, sep='\t', header=0)
    J1 = float(df.loc[df['directory index']==NumSearch+2,'QoI'])
    reduction = abs(J1-J0)/abs(J0+eps)

    df = pd.read_csv(CGLog, sep='\t', header=0)
    df.at[df.index[-1],'after linmin'] = J1
    df.at[df.index[-1],'reduction'] = reduction
    df.to_csv(CGLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)

    if (reduction<=tol):
        print ('FRPRMN - after linmin: conjugate-gradient optimization is finished.')
        commandFile = open(commandFilename,'w')
        commandFile.write('break\n')
        commandFile.close()
        commandFile = open(decisionMakerCommandFilename,'w')
        commandFile.write('break\n')
        commandFile.close()
        return 0

    #copy line minimization log
    import os
    numFiles = len(os.listdir('./linminLog/'))
    subprocess.check_call('cp '+lineMinLog+' linminLog/'+prefix+'.line_minimization.'+str(numFiles)+'.txt', shell=True)

    commandFile = open(commandFilename,'w')
    commandFile.write('cp '+str(NumSearch+2)+'/'+forwardFilename+' ./ \n')
    for k in range(NumCGFile):
        commandFile.write('cp '+str(NumSearch+2)+'/'+controlForcingFilenames[k]+' ./ \n')
        commandFile.write('mv '+gradientFilenames[k]+' '+'previous.'+gradientFilenames[k]+'\n')
        commandFile.write('mv '+CGFilenames[k]+' '+'previous.'+CGFilenames[k]+'\n')
    if(zeroBaseline):
        commandFile.write(magudiSetOptionCommand)
        commandFile.write('setOption "controller_switch" true\n')
    commandFile.write('srun -n '+str(NumProcs)+' ./forward\n')
    commandFile.write('srun -n '+str(NumProcs)+' ./adjoint\n')
    commandFile.close()
    commandFile = open(decisionMakerCommandFilename,'w')
    commandFile.write('python '+decisionMaker+' 1 \n')
    commandFile.close()

    print ('FRPRMN - after linmin: postprocessing is finished. Run new forward/adjoint simulations.')
    return 1
