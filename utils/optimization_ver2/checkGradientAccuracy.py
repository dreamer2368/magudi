import numpy as np
import subprocess
from base import *

command = forwardRunCommand()
subprocess.check_call(command, shell=True)

QoI0, subJ0 = QoI()
print "QoI0: ", QoI0

command = adjointRunCommand()
command += gatherControlForcingGradientCommand()
command += innerProductCommand(globalGradFiles,globalGradFiles,ggFiles)
subprocess.check_call(command, shell=True)

Grad0, subGrad0 = readInnerProduct(ggFiles)
print "Grad0: ", Grad0

Nk = 20
Ak = 10.0**(0.0-0.25*np.array(range(Nk)))
QoIk = np.zeros((Nk,),dtype=np.double)
Gradk = np.zeros((Nk,),dtype=np.double)
ek = np.zeros((Nk,),dtype=np.double)

stepFiles = globalControlSpaceFiles.copy()
for k, file in enumerate(stepFiles):
    stepFiles[k] = 'step/' + file

fId = open(prefix+'.gradient_accuracy.txt','w')
for k in range(Nk):
    command = zaxpyCommand(stepFiles,Ak[k],globalGradFiles,globalControlSpaceFiles)
    command += distributeCommand('step')
    command += forwardRunCommand('step')
    subprocess.check_call(command, shell=True)
    QoIk[k], subJ = QoI('step')
    Gradk[k] = (QoIk[k]-QoI0)/Ak[k]
    ek[k] = abs( (Gradk[k]-Grad)/Grad )
    print ("{:.16E}".format(Ak[k]), "{:.16E}".format(QoIk[k]), "{:.16E}".format(Gradk[k]), "{:.16E}".format(ek[k]))

    fId.write('%.16E\t%.16E\t%.16E\t%.16E' % (Ak[k], QoIk[k], Gradk[k], ek[k]))

fId.close()
