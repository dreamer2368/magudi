from base import *

def parabolic_interp(stepBracket, JBracket):
    a,b,c = stepBracket
    fa,fb,fc = JBracket

    x = b - 0.5 * ( ((b-a)**2)*(fb-fc) - ((b-c)**2)*(fb-fa) )                                           \
                / ( (b-a)*(fb-fc) - (b-c)*(fb-fa) )
    return x

def nextLinmin(zeroBaseline=True, initial=True, stop=False):
    import pandas as pd
    import subprocess

    if (initial):
        df = pd.read_csv(lineMinLog, sep='\t', header=0)
    else:
        df = collectQoI(lineMinLog)

        steps = np.array(df['step'][df['directory index']!='0'])
        Js = np.array(df['QoI'][df['directory index']!='0'])
        dirIdx = np.array(df['directory index'][df['directory index']!='0'])
        N = steps.size
        if (N!=4):
            raise RuntimeError('Cannot find values corresponding to a, b, c, and x!')
        minIdx = max(np.argmin(Js),1)

        temp = ['a','b','c']
        for k in range(3):
            switchDirectory(dirIdx[minIdx-1+k],temp[k],df)
            dirIdx = np.array(df['directory index'][df['directory index']!='0'])
        df.loc[df['directory index']=='x','directory index'] = '0'
        df.to_csv(lineMinLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)

    if (stop):
        print (df[df['directory index']!='0'])
        print ('LINMIN: line minimization is stopped.')

        commandFile = open(globalCommandFile,'w')
        commandFile.close()
        commandFile = open(decisionMakerCommandFile,'w')
        command = 'python3 '+decisionMaker+' 5'
        if(zeroBaseline):
            command += ' -zero_baseline'
        commandFile.write(command+'\n')
        commandFile.write('exit 0\n')
        commandFile.close()
        return 0

    a = np.array(df['step'][df['directory index']=='a'])[0]
    b = np.array(df['step'][df['directory index']=='b'])[0]
    c = np.array(df['step'][df['directory index']=='c'])[0]
    Ja = np.array(df['QoI'][df['directory index']=='a'])[0]
    Jb = np.array(df['QoI'][df['directory index']=='b'])[0]
    Jc = np.array(df['QoI'][df['directory index']=='c'])[0]

    if ( (c-a) < b * linminTol + eps ):
        print ('steps: %.16E %.16E %.16E'%(a,b,c))
        print ('QoIs: %.16E %.16E %.16E'%(Ja,Jb,Jc))
        print ('LINMIN: line minimization is finished.')

        commandFile = open(globalCommandFile,'w')
        commandFile.close()
        commandFile = open(decisionMakerCommandFile,'w')
        command = 'python3 '+decisionMaker+' 5'
        if(zeroBaseline):
            command += ' -zero_baseline'
        commandFile.write(command+'\n')
        commandFile.close()
        return 0

    xs = parabolic_interp([a,b,c],[Ja,Jb,Jc])
    Cr = 1.0 - 1.0/golden_ratio
    if( (xs>c) or (xs<a) or (abs(np.log10((c-b)/(b-a)))>1.0) ):
        if( b>0.5*(a+c) ):
            xs = b - Cr * (b-a)
        else:
            xs = b + Cr * (c-b)
    data = {'step':[xs], 'QoI':[np.nan], 'directory index':['x']}
    new_df = pd.DataFrame(data)
    df = df.append(new_df, ignore_index=True)

    commandString = ''
    temp = globalControlSpaceFiles.copy()
    if (zeroBaseline):
        for j in range(NcontrolRegion):
            temp[j] = ''
    target = globalControlSpaceFiles.copy()
    for k, file in enumerate(target):
        target[k] = 'x/' + file
        commandString += 'rm %s \n'%target[k]
    commandString += zaxpyCommand(target, -xs, globalConjugateGradientFiles, temp)
    commandString += '\n'
    commandString += forwardRunCommand('x')
    commandString += bashCheckResultCommand('intermediate forward runs')
    fID = open(globalCommandFile,'w')
    fID.write(commandString)
    fID.close()

    commandFile = open(decisionMakerCommandFile,'w')
    command = 'python3 '+decisionMaker+' 4'
    if(zeroBaseline):
        command += ' -zero_baseline'
    commandFile.write(command+'\n')
    commandFile.close()

    df.to_csv(lineMinLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)
    print (df[df['directory index']!='0'])
    print ('LINMIN: next linmin evaluation is prepared-')
    return 1
