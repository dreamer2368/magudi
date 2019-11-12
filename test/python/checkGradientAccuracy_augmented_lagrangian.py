import numpy as np
import subprocess
from base import *

commandFile = 'test.sh'

#fID = open(commandFile,'w')
#command = forwardRunCommand('.',True)
#fID.write(command)
#fID.close()
#subprocess.check_call('bash '+commandFile, shell=True)

QoI0, subJ0 = QoI()
print ("QoI0: ", QoI0, subJ0)

#from os import path
#for file in globalGradFiles:
#    if(path.exists(file)):
#        subprocess.check_call('rm '+file,shell=True)
#
#fID = open(commandFile,'w')
#command = adjointRunCommand()
#command += gatherControlForcingGradientCommand()
#command += innerProductCommand(globalGradFiles,globalGradFiles,ggFiles)
#fID.write(command)
#fID.close()
#subprocess.check_call('bash '+commandFile, shell=True)

Grad0, subGrad0 = readInnerProduct(ggFiles)
print ("Grad0: ", Grad0, subGrad0)

Nk = 20
Ak = 10.0**(2.0-0.25*np.array(range(Nk)))
QoIk = np.zeros((Nk,),dtype=np.double)
Gradk = np.zeros((Nk,),dtype=np.double)
ek = np.zeros((Nk,),dtype=np.double)

baselineFiles = globalControlSpaceFiles.copy()
baselineFiles[0] = ''
stepFiles = globalControlSpaceFiles.copy()
for k, file in enumerate(stepFiles):
    stepFiles[k] = 'step/' + file

fId = open('test.gradient_accuracy.txt','w')
for k in range(Nk):
    command = zaxpyCommand(stepFiles,Ak[k],globalGradFiles,baselineFiles)
    command += forwardRunCommand('step',False)
    fID = open(commandFile,'w')
    fID.write(command)
    fID.close()
    subprocess.check_call('bash '+commandFile, shell=True)
    QoIk[k], subJ = QoI('step')
    Gradk[k] = (QoIk[k]-QoI0)/Ak[k]
    ek[k] = abs( (Gradk[k]-Grad0)/Grad0 )
    print ("{:.16E}".format(Ak[k]), "{:.16E}".format(QoIk[k]), "{:.16E}".format(Gradk[k]), "{:.16E}".format(ek[k]))

    fId.write('%.16E\t%.16E\t%.16E\t%.16E\n' % (Ak[k], QoIk[k], Gradk[k], ek[k]))

fId.close()
