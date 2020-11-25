import numpy as np
import subprocess
from base import *

command = forwardRunCommand('x0',True)
fID = open('test-forward.sh','w')
fID.write(command)
fID.close()
subprocess.check_call('bash test-forward.sh', shell=True)

command = adjointRunCommand()
command += gatherControlForcingGradientCommand()
command += innerProductCommand(globalGradFiles,globalGradFiles,ggFiles)
fID = open('test-adjoint.sh','w')
fID.write(command)
fID.close()
subprocess.check_call('bash test-adjoint.sh', shell=True)

QoI0, subJ0 = QoI('x0')
print ("QoI0: ", QoI0)
print ("subJ: ", subJ0)

Grad0, subGrad0 = readInnerProduct(ggFiles)
print ("Grad0: ", Grad0)
print ("subGrad: ", subGrad0)

Nk = 20
Ak = 10.0**(-2.0-0.25*np.array(range(Nk)))
QoIk = np.zeros((Nk,),dtype=np.double)
Gradk = np.zeros((Nk,),dtype=np.double)
ek = np.zeros((Nk,),dtype=np.double)

temp = ['x0/'+file for file in globalControlSpaceFiles]
for k in range(NcontrolRegion):
    temp[k] = ''

stepFiles = globalControlSpaceFiles.copy()
for k, file in enumerate(stepFiles):
    stepFiles[k] = 'x/' + file

for k in range(Nk):
    amp = Ak[k]
    command = zaxpyCommand(stepFiles,amp,globalGradFiles,temp)
    command += forwardRunCommand('x')
    command += scriptor.purgeDirectoryCommand('x')
    fID = open('test-step.sh','w')
    fID.write(command)
    fID.close()
    subprocess.check_call('bash test-step.sh', shell=True)
    QoIk[k], subJ = QoI('x')
    Gradk[k] = (QoIk[k]-QoI0)/Ak[k]
    ek[k] = abs( (Gradk[k]-Grad0)/Grad0 )
    print ("{:.16E}".format(Ak[k]), "{:.16E}".format(QoIk[k]), "{:.16E}".format(Gradk[k]), "{:.16E}".format(ek[k]))

    fId = open(globalPrefix+'.gradient_accuracy.txt','a+')
    fId.write('%.16E\t%.16E\t%.16E\t%.16E\n' % (Ak[k], QoIk[k], Gradk[k], ek[k]))
    fId.close()
