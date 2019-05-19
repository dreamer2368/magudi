NumProcs = 360
NumSearch = 10
prefix = 'MultiblockJet'
lineMinLog = prefix+'.line_minimization.txt'
CGLog = prefix+'.conjugate_gradient.txt'
dggFilename = prefix+'.dgg.txt'

forwardFilename = prefix+'.forward_run.txt'
adjointFilename = prefix+'.adjoint_run.txt'
controlForcingFilenames = [prefix+'.control_forcing_controlRegion.E.dat',                       \
                            prefix+'.control_forcing_controlRegion.W.dat',                      \
                            prefix+'.control_forcing_controlRegion.N.dat',                      \
                            prefix+'.control_forcing_controlRegion.S.dat']
gradientFilenames = [prefix+'.gradient_controlRegion.E.dat',                                    \
                     prefix+'.gradient_controlRegion.W.dat',                                    \
                     prefix+'.gradient_controlRegion.N.dat',                                    \
                     prefix+'.gradient_controlRegion.S.dat']
CGFilenames = [prefix+'.conjugate_gradient_controlRegion.E.dat',                                \
                prefix+'.conjugate_gradient_controlRegion.W.dat',                               \
                prefix+'.conjugate_gradient_controlRegion.N.dat',                               \
                prefix+'.conjugate_gradient_controlRegion.S.dat']
normFilenames = [prefix+'.norm_controlRegion.E.dat',                                            \
                 prefix+'.norm_controlRegion.W.dat',                                            \
                 prefix+'.norm_controlRegion.N.dat',                                            \
                 prefix+'.norm_controlRegion.S.dat']

decisionMaker = 'optimization.py'
commandFilename = prefix+'.command.sh'
decisionMakerCommandFilename = prefix+'.command.python.ready.sh'

initial_step = 0.1
golden_ratio = 1.618034
tol, eps = 1.0e-7,  1.0e-7

magudiSetOptionCommand = 'function setOption() {\n'                                         \
                         '    if grep -q "$1" magudi.inp\n'                                 \
                         '    then\n'                                                       \
                         '	sed -i "s/^.*$1.*$/$1 = $2/g" magudi.inp\n'                     \
                         '    else\n'                                                       \
                         '	echo "$1 = $2" >> magudi.inp\n'                                 \
                         '    fi\n'                                                         \
                         '}\n'

import numpy as np
import subprocess
import pandas as pd

def bashCheckResultCommand(procedureName):
    commandString = 'if [ $? -ne 0 ]; then\n'                                               \
                    '   echo "'+procedureName+' failed."\n'                                 \
                    '   exit -1\n'                                                          \
                    'else\n'                                                                \
                    '   echo "'+procedureName+' succeeded."\n'                              \
                    '   exit 0\n'                                                           \
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
