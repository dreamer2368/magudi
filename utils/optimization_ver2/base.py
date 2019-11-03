from constants import *
from filenames import *

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
    
def QoI(rootDirectory = None):
    if (rootDirectory is None):
        rdir = '.'
    else:
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
    for j in range(NcontrolRegion):
        command = './zxdoty %s %s %s %s \n'                                                        \
                % (outputFiles_[j], xFiles[j], yFiles[j], globalNormFiles[j])
        commandString += command

    for k in range(NcontrolRegion,NcontrolRegion+Nsplit-1):
        index = k+1 - NcontrolRegion
        inputFile = '%s/%s'%(directories[index],inputFiles[index])
        commandString += './spatial_inner_product %s %s --input %s --output %s \n'                 \
                            % (xFiles[k], yFiles[k], inputFile, outputFiles_[k])

    return commandString

def readInnerProduct(outputFiles_=None):
    if (outputFiles_ is None):
        outputFiles_ = innerProductFiles
        
    scalar = 0.0
    subScalar = np.zeros(NcontrolRegion+Nsplit-1)
    for k in range(NcontrolRegion+Nsplit-1):
        subScalar[k] = readScalar(outputFiles_[k])
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
    for j in range(NcontrolRegion):
        commandString += './zaxpy ' + zFiles[j] + ' ' + "{:.16E}".format(a) + ' ' + xFiles[j]
        if (yFiles is not None):
            commandString += ' '+yFiles[j]+'\n'
        else:
            commandString += '\n'
        commandString += bashCheckResultCommand('zaxpy-control-forcing-%d'%j)
            
    w2 = initialConditionControllability
    for k in range(NcontrolRegion,NcontrolRegion+Nsplit-1):
        index = k+1 - NcontrolRegion
        inputFile = '%s/%s'%(directories[index],inputFiles[index])
        commandString += './qfile_zaxpy ' + zFiles[k] + ' ' + "{:.16E}".format(a*w2) + ' ' + xFiles[k]
        if (yFiles is not None):
            commandString += ' '+yFiles[k]+'\n'
        else:
            commandString += '\n'
        commandString += bashCheckResultCommand('zaxpy-initial-condition-%d'%index)
            
    return commandString

def distributeCommand(rootDirectory):
    rdir = rootDirectory
    
    commandString = ''
    for k in range(Nsplit):
        kOffset = NtimestepOffset * k
        for j in range(NcontrolRegion):
            sliceControlForcingFile = '%s/%s/%s%s'%(rdir,directories[k],prefixes[k],controlForcingFiles[j])
            globalControlForcingFile = '%s/%s'%(rdir,globalControlSpaceFiles[j])
            commandString += './slice_control_forcing ' + sliceControlForcingFile + ' ' + globalControlForcingFile + ' ' + \
                             str(int(totalTimestep)) + ' ' + str(int(kOffset)) + ' ' + str(int(Nts)) + '\n'
        
    for k in range(Nsplit):
        commandString += 'cp ' + '%s/%s'%(rdir,icFiles[k]) + ' ' + '%s/%s/%s'%(rdir,directories[k],icFiles[k]) + '\n'
        
    return commandString

def gatherControlForcingGradientCommand():
    commandString = ''
    for k in range(Nsplit):
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

def forwardRunCommand(rootDirectory=None):
    if (rootDirectory is None):
        rdir = '.'
    else:
        rdir = rootDirectory
    
    commandString = ''
    commandString += distributeCommand(rdir)
    
    for k in range(Nsplit):
        commandString += 'cd %s/%s \n'%(rdir,directories[k])
        commandString += './forward --input %s \n' % inputFiles[k]
        if (rdir=='.'):
            commandString += 'cd .. \n'
        else:
            commandString += 'cd ../../ \n'
            
    for k in range(Nsplit-1):
        kOffset = startTimestep + NtimestepOffset * (k+1)
        matchingFile = '%s/%s/%s'%(rdir,directories[k],matchingForwardFiles[k])
        icFile = '%s/%s/%s'%(rdir,directories[k+1],icFiles[k+1])
        diffFile = '%s/%s'%(rdir,diffFiles[k])
        inputFile = '%s/%s/%s'%(rdir,directories[k],inputFiles[k])
        commandString += './qfile_zaxpy %s %.16E %s %s --input %s \n'                       \
                            % (diffFile,-1.0,icFile,matchingFile,inputFile)
        
    for k in range(Nsplit-1):
        diffFile = '%s/%s'%(rdir,diffFiles[k])
        inputFile = '%s/%s/%s'%(rdir,directories[k],inputFiles[k])
        diffOutputFile = '%s/%s'%(rdir,diffOutputFiles[k])
        commandString += './spatial_inner_product %s %s --output %s --input %s \n'          \
                            % (diffFile, diffFile, diffOutputFile, inputFile)
    return commandString

def adjointRunCommand(rootDirectory=None):
    if (rootDirectory is None):
        rdir = '.'
    else:
        rdir = rootDirectory

    commandString = ''
    for k in range(Nsplit-1):
        commandString += 'cd %s/%s \n' % (rdir,directories[k])
        commandString += setOptionCommand(inputFiles[k])
        commandString += 'setOption "enable_adjoint_restart" "true"\n'
        commandString += 'setOption "number_of_timesteps" %d \n' % (Nts-NtimestepOffset)
        commandString += 'setOption "adjoint_restart\/accumulated_timesteps" 0 \n'
        commandString += 'setOption "adjoint_restart\/intermediate_end_timestep" %d \n' % (NtimestepOffset)
        commandString += './adjoint --input %s \n' % inputFiles[k]
        if (rdir=='.'):
            commandString += 'cd .. \n'
        else:
            commandString += 'cd ../../ \n'
        commandString += '\n'
            
    k = Nsplit - 1
    commandString += 'cd %s/%s \n' % (rdir,directories[k])
    commandString += setOptionCommand(inputFiles[k])
    commandString += 'setOption "enable_adjoint_restart" "false"\n'
    commandString += 'setOption "number_of_timesteps" %d \n' % (Nts)
    commandString += './adjoint --input %s \n' % inputFiles[k]
    if (rdir=='.'):
        commandString += 'cd .. \n'
    else:
        commandString += 'cd ../../ \n'
    commandString += '\n'
    
    for k in range(Nsplit-1):
        matchingAdjointFile = '%s/%s'%(directories[k],matchingAdjointFiles[k])
        inputFile = '%s/%s/%s'%(rdir,directories[k],inputFiles[k])
        diffFile = '%s/%s'%(rdir,diffFiles[k])
        commandString += './qfile_zaxpy %s %.16E %s %s --input %s \n'                              \
        % (matchingAdjointFile, matchingConditionWeight, diffFile, matchingAdjointFile, inputFile)
        commandString += '\n'
    
    for k in range(Nsplit-1):
        commandString += 'cd %s/%s \n' % (rdir,directories[k])
        commandString += setOptionCommand(inputFiles[k])
        commandString += '\n'
        commandString += 'setOption "number_of_timesteps" %d \n' % (NtimestepOffset)
        commandString += 'setOption "adjoint_restart\/accumulated_timesteps" %d \n' % (Nts-NtimestepOffset)
        commandString += 'setOption "adjoint_restart\/intermediate_end_timestep" 0 \n'
        commandString += './adjoint --input %s \n' % inputFiles[k]
        commandString += 'setOption "number_of_timesteps" %d \n' % (Nts)
        commandString += 'setOption "enable_adjoint_restart" "false"\n'
        if (rdir=='.'):
            commandString += 'cd .. \n'
        else:
            commandString += 'cd ../../ \n'
        commandString += '\n'
    
    for k in range(1,Nsplit):
        icAdjointFile = '%s/%s/%s' % (rdir,directories[k],icAdjointFiles[k])
        diffFile = '%s/%s.diff.q'%(rdir,prefixes[k-1])
        icGradFile = '%s/%s' % (rdir,icGradientFiles[k])
        inputFile = '%s/%s/%s'%(rdir,directories[k],inputFiles[k])
        commandString += './qfile_zaxpy %s %.16E %s %s --input %s \n'                             \
        % (icGradFile, -matchingConditionWeight, diffFile, icAdjointFile, inputFile)
        commandString += '\n'

    return commandString

def dggCommand():
    commandString = ''
    for j in range(NcontrolRegion):
        commandString += './zwxmwy %s %s %s previous.%s %s \n'                                            \
        %(dggFiles[j], globalGradFiles[j], globalGradFiles[j], globalGradFiles[j], globalNormFiles[j])
    commandString += '\n'
    for k in range(1,Nsplit):
        diffFile = '%s.diff.adjoint.q'%(prefixes[k])
        icGradFile = '%s' % (icGradientFiles[k])
        inputFile = '%s/%s'%(directories[k],inputFiles[k])
        commandString += './qfile_zaxpy %s %.16E previous.%s %s --input %s \n'                             \
                            % (diffFile, -1.0, icGradFile, icGradFile, inputFile)
    commandString += '\n'
    for k in range(1,Nsplit):
        idx = NcontrolRegion-1 + k
        diffFile = '%s.diff.adjoint.q'%(prefixes[k])
        inputFile = '%s/%s'%(directories[k],inputFiles[k])
        commandString += './spatial_inner_product %s %s --output %s --input %s \n'          \
                            % (diffFile, icGradientFiles[k], dggFiles[idx], inputFile)
    return commandString
        