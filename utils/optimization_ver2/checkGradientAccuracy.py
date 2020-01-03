import numpy as np
import subprocess
from base import *

command = forwardRunCommand('.',True)
fID = open('test-forward.sh','w')
fID.write(command)
fID.close()
subprocess.check_call('bash test-forward.sh', shell=True)

QoI0, subJ0 = QoI()
print ("QoI0: ", QoI0)
print ("subJ: ", subJ0)

command = adjointRunCommand()
command += gatherControlForcingGradientCommand()
command += innerProductCommand(globalGradFiles,globalGradFiles,ggFiles)
fID = open('test-adjoint.sh','w')
fID.write(command)
fID.close()
subprocess.check_call('bash test-adjoint.sh', shell=True)

Grad0, subGrad0 = readInnerProduct(ggFiles)
print ("Grad0: ", Grad0)
print ("subGrad: ", subGrad0)

Nk = 20
Ak = 10.0**(3.0-0.25*np.array(range(Nk)))
QoIk = np.zeros((Nk,),dtype=np.double)
Gradk = np.zeros((Nk,),dtype=np.double)
ek = np.zeros((Nk,),dtype=np.double)

temp = globalControlSpaceFiles.copy()
for k in range(NcontrolRegion):
    temp[k] = ''

stepFiles = globalControlSpaceFiles.copy()
for k, file in enumerate(stepFiles):
    stepFiles[k] = 'step/' + file

fId = open(globalPrefix+'.gradient_accuracy.txt','w')
for k in range(Nk):
#for k in range(1):
    amp = Ak[k]
    #amp = 0.0
    command = zaxpyCommand(stepFiles,amp,globalGradFiles,temp)
    command += forwardRunCommand('step')
    fID = open('test-step.sh','w')
    fID.write(command)
    fID.close()
    subprocess.check_call('bash test-step.sh', shell=True)
    QoIk[k], subJ = QoI('step')
    Gradk[k] = (QoIk[k]-QoI0)/Ak[k]
    ek[k] = abs( (Gradk[k]-Grad0)/Grad0 )
    print ("{:.16E}".format(Ak[k]), "{:.16E}".format(QoIk[k]), "{:.16E}".format(Gradk[k]), "{:.16E}".format(ek[k]))

    fId.write('%.16E\t%.16E\t%.16E\t%.16E' % (Ak[k], QoIk[k], Gradk[k], ek[k]))

fId.close()
