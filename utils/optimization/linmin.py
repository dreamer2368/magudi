from base import *

def parabolic_interp(stepBracket, JBracket):
    a,b,c = stepBracket
    fa,fb,fc = JBracket

    x = b - 0.5 * ( ((b-a)**2)*(fb-fc) - ((b-c)**2)*(fb-fa) )                                           \
                / ( (b-a)*(fb-fc) - (b-c)*(fb-fa) )
    return x

def nextLinmin(forwardFilename, CGFilenames, controlForcingFilenames, zeroBaseline=True, initial=True, stop=False):
    import pandas as pd
    import subprocess
    from base import readScalar
    # from IPython.display import display, HTML
    NumCGFile = len(CGFilenames)
    if (len(controlForcingFilenames)!=NumCGFile):
        print ('LINMIN : numbers of control space files do not match!')
        return

    if (initial):
        df = pd.read_csv(lineMinLog, sep='\t', header=0)
    else:
        df = collectQoIs(lineMinLog, forwardFilename, NumSearch)

        steps = np.array(df['step'][df['directory index']>0])
        Js = np.array(df['QoI'][df['directory index']>0])
        dirIdx = np.array(df['directory index'][df['directory index']>0])
        N = steps.size
        minIdx = max(np.argmin(Js),1)

        for k in range(3):
            switchDirectory(dirIdx[minIdx-1+k],NumSearch+1+k,df)
            dirIdx = np.array(df['directory index'][df['directory index']>0])
#         display(df)
        df.to_csv(lineMinLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)

    if (stop):
        print (df[df['directory index']>0])
        print ('LINMIN: line minimization is stopped.')

        commandFile = open(commandFilename,'w')
        commandFile.close()
        commandFile = open(decisionMakerCommandFilename,'w')
        command = 'python3 '+decisionMaker+' 5'
        if(zeroBaseline):
            command += ' -zero_baseline'
        commandFile.write(command+'\n')
        commandFile.write('exit 0\n')
        commandFile.close()
        return 0

    stepBracket = np.array(df['step'][df['directory index']>NumSearch])
    JBracket = np.array(df['QoI'][df['directory index']>NumSearch])
    if (stepBracket.size!=3):
        print ('LINMIN: the minimum bracket is lost or not found yet!')
        return -1
    bracketSize = (stepBracket[2]-stepBracket[0])/2.0/(stepBracket[1]+eps)
    if (bracketSize<=1.0e-1):
        print (stepBracket)
        print (JBracket)
        print ('LINMIN: line minimization is finished.')

        commandFile = open(commandFilename,'w')
        commandFile.close()
        commandFile = open(decisionMakerCommandFilename,'w')
        command = 'python3 '+decisionMaker+' 5'
        if(zeroBaseline):
            command += ' -zero_baseline'
        commandFile.write(command+'\n')
#        commandFile.write('exit 1\n')
        commandFile.close()
        return 0

    xs = parabolic_interp(stepBracket,JBracket)
    if (np.remainder(NumSearch,2)==1):
        idx = int(NumSearch/2)
        nLeft, nRight = int(NumSearch/2), int(NumSearch/2)
    elif ( xs < 0.5*(stepBracket[0]+stepBracket[2]) ):
        idx = int(NumSearch/2) - 1
        nLeft, nRight = int(NumSearch/2)-1, int(NumSearch/2)
    else:
        idx = int(NumSearch/2)
        nLeft, nRight = int(NumSearch/2), int(NumSearch/2)-1

    steps, Js = np.zeros(NumSearch), np.zeros(NumSearch)
    steps[idx] = xs
    if (stepBracket[0]>0.0):
        dx = xs * ( (stepBracket[0]/xs)**( 1./(nLeft+1)*np.arange(1,nLeft+1) ) )
        dx -= stepBracket[0]
        steps[:nLeft] = xs - dx
    else:
        dx = xs * ( golden_ratio**( -nLeft-1 + np.arange(1,nLeft+1) ) )
        steps[:nLeft] = dx
    steps[-nRight:] = xs * ( (stepBracket[2]/xs)**( 1./(nRight+1)*np.arange(1,nRight+1) ) )

    df.loc[df['directory index']<=NumSearch,'directory index'] = 0
    data = {'step':steps, 'QoI':np.ones(NumSearch)*np.nan, 'directory index':list(range(1,NumSearch+1))}
    new_df = pd.DataFrame(data)
    df = df.append(new_df, ignore_index=True)

    commandFile = open(commandFilename,'w')
    for k in range(NumSearch):
        command = 'rm '+str(k+1)+'/'+prefix+'.*.dat'
        commandFile.write(command+'\n')
        for i in range(NumCGFile):
            command = 'srun -n '+str(NumProcs)+' ./zaxpy '                                                 \
                        +str(k+1)+'/'+controlForcingFilenames[i]+' '                                    \
                        +"{:.16E}".format(-steps[k])+' '+CGFilenames[i]
            if (not zeroBaseline):
                command += ' '+controlForcingFilenames[i]
            commandFile.write(command+'\n')
            commandFile.write(bashCheckResultCommand('zaxpy-'+str(k)+'-'+str(i)))
    commandFile.write('bash intermediate_forward_runs.sh\n')
    commandFile.write(bashCheckResultCommand('intermediate forward runs'))
    commandFile.close()
    commandFile = open(decisionMakerCommandFilename,'w')
    command = 'python3 '+decisionMaker+' 4'
    if(zeroBaseline):
        command += ' -zero_baseline'
    commandFile.write(command+'\n')
    commandFile.close()

    df.to_csv(lineMinLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)
    print ('submitted all zaxpy works. Wait until they finish.')
    print (df[df['directory index']>0])
    print ('LINMIN: next linmin evaluation is prepared-')
    return 1
