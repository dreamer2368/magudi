NumSearch = 10
Nsplit = 3                                               # number of split time segments
Nts = 6                                                  # number of timesteps of each time segment
NtimestepOffset = 4                                      # timestep offset between time segments
startTimestep = 12                                       # initial timestep of the first time segment
totalTimestep = (Nsplit-1) * NtimestepOffset + Nts       # number of timesteps for the entire time span

matchingConditionWeight = 1.0e-7                         # weight for matching condition penalty
initialConditionControllability = 1.0e-9                 # weight for derivative with respect to initial conditions

NcontrolRegion = 1                                       # number of control region

initial_step = 0.1
golden_ratio = 1.618034
tol, eps = 1.0e-7,  1.0e-7

import numpy as np
import subprocess
import pandas as pd

globalPrefix = 'AcousticMonopole'
prefixes, inputFiles, outputFiles, commandFiles = [], [], [], []
icFiles, diffFiles, diffOutputFiles             = [], [], []
matchingForwardFiles, matchingAdjointFiles      = [], []
icAdjointFiles, icGradientFiles, directories    = [], [], []
for k in range(Nsplit):
    prefixes += ['%s-%01d'%(globalPrefix,k)]
    inputFiles += ['magudi-%01d.inp'%(k)]
    outputFiles += ['%s-%01d.forward_run.txt'%(globalPrefix,k)]
    commandFiles += ['%s-%01d.command.sh'%(globalPrefix,k)]
    icFiles += ['%s-%01d.ic.q'%(globalPrefix,k)]
    
    kOffset = startTimestep + NtimestepOffset * (k+1)
    diffFiles += ['%s-%01d.diff.q'%(globalPrefix,k)]                             #use only 0, ... , Nsplit-2
    diffOutputFiles += ['%s-%01d.diff.txt'%(globalPrefix,k)]                     #use only 0, ... , Nsplit-2
    matchingForwardFiles += ['%s-%01d-%08d.q'%(globalPrefix,k,kOffset)]          #use only 0, ... , Nsplit-2
    matchingAdjointFiles += ['%s-%01d-%08d.adjoint.q'%(globalPrefix,k,kOffset)]  #use only 0, ... , Nsplit-2
    
    kOffset = startTimestep + NtimestepOffset * k
    icAdjointFiles += ['%s-%01d-%08d.adjoint.q'%(globalPrefix,k,kOffset)]        #use only 1, ... , Nsplit-1
    icGradientFiles += ['%s-%01d.ic.adjoint.q'%(globalPrefix,k)]                 #use only 1, ... , Nsplit-1
    directories += ['%01d'%k]
    
diffFiles[-1] = ''
diffOutputFiles[-1] = ''
matchingForwardFiles[-1] = ''
matchingAdjointFiles[-1] = ''
icAdjointFiles[0] = ''
icGradientFiles[0] = ''

controlRegions = ['controlRegion']
gradientFiles, normFiles, conjugateGradientFiles, controlForcingFiles = [], [], [], []
for k in range(NcontrolRegion):
    gradientFiles += ['.gradient_%s.dat'%controlRegions[k]]
    normFiles += ['.norm_%s.dat'%controlRegions[k]]
    conjugateGradientFiles += ['.conjugate_gradient_%s.dat'%controlRegions[k]]
    controlForcingFiles += ['.control_forcing_%s.dat'%controlRegions[k]]
    
globalGradFiles, globalNormFiles, globalConjugateGradientFiles, globalControlForcingFiles = [], [], [], []
globalSensitivityFiles = []
for k in range(NcontrolRegion):
    globalGradFiles += [globalPrefix + gradientFiles[k]]
    globalNormFiles += [globalPrefix + normFiles[k]]
    globalConjugateGradientFiles += [globalPrefix + conjugateGradientFiles[k]]
    globalControlForcingFiles += [globalPrefix + controlForcingFiles[k]]
    globalSensitivityFiles += ['%s.adjoint_run_%s.txt'%(globalPrefix,controlRegions[k])]
    
innerProductFiles = []
for j in range(NcontrolRegion):
    innerProductFiles += ['%s.inner_product_%s.txt'%(globalPrefix,controlRegions[j])]
for k in range(1,Nsplit):
    innerProductFiles += ['%s-%01d.inner_product.txt'%(globalPrefix,k)]
    
globalCommandFile = globalPrefix + '.command.sh'

# lineMinLog = prefix+'.line_minimization.txt'
# CGLog = prefix+'.conjugate_gradient.txt'
# dggFilename = prefix+'.dgg.txt'

# forwardFilename = prefix+'.forward_run.txt'
# adjointFilename = prefix+'.adjoint_run.txt'
# controlForcingFilenames = [prefix+'.control_forcing_controlRegion.E.dat',                       \
#                             prefix+'.control_forcing_controlRegion.W.dat',                      \
#                             prefix+'.control_forcing_controlRegion.N.dat',                      \
#                             prefix+'.control_forcing_controlRegion.S.dat']
# gradientFilenames = [prefix+'.gradient_controlRegion.E.dat',                                    \
#                      prefix+'.gradient_controlRegion.W.dat',                                    \
#                      prefix+'.gradient_controlRegion.N.dat',                                    \
#                      prefix+'.gradient_controlRegion.S.dat']
# CGFilenames = [prefix+'.conjugate_gradient_controlRegion.E.dat',                                \
#                 prefix+'.conjugate_gradient_controlRegion.W.dat',                               \
#                 prefix+'.conjugate_gradient_controlRegion.N.dat',                               \
#                 prefix+'.conjugate_gradient_controlRegion.S.dat']
# normFilenames = [prefix+'.norm_controlRegion.E.dat',                                            \
#                  prefix+'.norm_controlRegion.W.dat',                                            \
#                  prefix+'.norm_controlRegion.N.dat',                                            \
#                  prefix+'.norm_controlRegion.S.dat']

decisionMaker = 'optimization.py'
# commandFilename = prefix+'.command.sh'
# decisionMakerCommandFilename = prefix+'.command.python.ready.sh'

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

def bashCheckResultCommand(procedureName):
    commandString = 'if [ $? -ne 0 ]; then\n'                                               \
                    '   echo "'+procedureName+' failed."\n'                                 \
                    '   exit -1\n'                                                          \
                    'else\n'                                                                \
                    '   echo "'+procedureName+' succeeded."\n'                              \
                    'fi\n'
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
    
def QoI(rootDirectory):
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

def innerProductCommand(xFiles,yFiles):
    Nx, Ny = len(xFiles), len(yFiles)
    if ((Nx!=Nsplit+NcontrolRegion-1) or (Ny!=Nsplit+NcontrolRegion-1)):
        raise LookupError('Number of files must be equal to '+str(Nsplit+NcontrolRegion-1)+'!')
    
    commandString = ''
    for j in range(NcontrolRegion):
        command = './zxdoty ' + innerProductFiles[j] + ' ' +                                \
                  xFiles[j] + ' ' + yFiles[j] + ' ' + globalNormFiles[j] + '\n'
        commandString += 'echo ' + command
        commandString += command

    for k in range(NcontrolRegion,NcontrolRegion+Nsplit-1):
        index = k+1 - NcontrolRegion
        inputFile = '%s/%s'%(directories[index],inputFiles[index])
        commandString += './spatial_inner_product '+xFiles[k]+' '+yFiles[k] +                      \
                         ' --input ' + inputFile + ' --output ' + innerProductFiles[k] + '\n'

    return commandString

def readInnerProduct(outputFile=None):
    scalar = 0.0
    subScalar = np.zeros(NcontrolRegion+Nsplit-1)
    for k in range(NcontrolRegion+Nsplit-1):
        subScalar[k] = readScalar(innerProductFiles[k])
        scalar += subScalar[k]
        
    if(outputFile is not None):
        fID = open(outputFile,'w')
        fID.write('{:.15E}'.format(scalar))
        fID.close()
        
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
    for j in range(NcontrolRegion):
        commandString += './zaxpy ' + zFiles[j] + ' ' + "{:.16E}".format(a) + ' ' + xFiles[j]
        if (yFiles is not None):
            commandString += ' '+yFiles[j]+'\n'
        else:
            commandString += '\n'
            
    w2 = initialConditionControllability
    for k in range(NcontrolRegion,NcontrolRegion+Nsplit-1):
        index = k+1 - NcontrolRegion
        inputFile = '%s/%s'%(directories[index],inputFiles[index])
        commandString += './qfile_zaxpy ' + zFiles[k] + ' ' + "{:.16E}".format(a*w2) + ' ' + xFiles[k]
        if (yFiles is not None):
            commandString += ' '+yFiles[k]+'\n'
        else:
            commandString += '\n'
            
    return commandString

def distributeCommand(rootDirectory):
    rdir = rootDirectory
    
    commandString = ''
    for k in range(Nsplit):
        kOffset = NtimestepOffset * k
        for j in range(NcontrolRegion):
            sliceControlForcingFile = '%s/%s/%s%s'%(rdir,directories[k],prefixes[k],controlForcingFiles[j])
            globalControlForcingFile = '%s/%s'%(rdir,globalControlForcingFiles[j])
            commandString += './slice_control_forcing ' + sliceControlForcingFile + ' ' + globalControlForcingFile + ' ' + \
                             str(int(totalTimestep)) + ' ' + str(int(kOffset)) + ' ' + str(int(Nts)) + '\n'
        commandString += 'cp ' + '%s/%s'%(rdir,icFiles[k]) + ' ' + '%s/%s/%s'%(rdir,directories[k],icFiles[k]) + '\n'
        
    return commandString

def gatherControlForcingGradientCommand():
    commandString = ''
    for k in range(0):
        kOffset = NtimestepOffset * k
        for j in range(NcontrolRegion):
            sliceGradFile = '%s/%s%s'%(directories[k],prefixes[k],gradientFiles[j])
            commandString += './paste_control_forcing ' + globalGradFiles[j] + ' ' + sliceGradFile + ' ' + \
                              str(int(totalTimestep)) + ' ' + str(int(kOffset)) + ' ' + str(int(Nts)) + '\n'
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

def collectQoIs(logFilename, forwardFilename, NumSearch):
    import pandas as pd

    df = pd.read_csv(logFilename, sep='\t', header=0)

    for k in range(1,NumSearch+1):
        directory = str(k)+'/'
        df.at[df['directory index']==k,'QoI'] = readScalar(directory+forwardFilename)
    df = df.sort_values(by='step',ascending=True)
    df = df.reset_index(drop=True)
    return df
