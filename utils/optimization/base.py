NumProcs = 36
NumSearch = 10
prefix = 'AcousticMonopole'
lineMinLog = prefix+'.line_minimization.txt'
CGLog = prefix+'.conjugate_gradient.txt'
dggFilename = prefix+'.dgg.txt'

forwardFilename = prefix+'.forward_run.txt'
adjointFilename = prefix+'.adjoint_run.txt'
controlForcingFilenames = [prefix+'.control_forcing_controlRegion.dat']
gradientFilenames = [prefix+'.gradient_controlRegion.dat']
CGFilenames = [prefix+'.conjugate_gradient_controlRegion.dat']
normFilenames = [prefix+'.norm_controlRegion.dat']

initial_step = 0.1
golden_ratio = 1.618034
tol, eps = 1.0e-7,  1.0e-7

import numpy as np
import subprocess
import pandas as pd

def readScalar(scalarFilename):
    fID = open(scalarFilename,'r')
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
