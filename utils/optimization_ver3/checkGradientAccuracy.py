import numpy as np
import subprocess
from inputs import InputParser
from optimizer import Optimizer

def checkFiniteDifferenceOrder(filename):
    result = np.loadtxt(filename)

    order = np.log(result[1:,-1] / result[:-1,-1]) / np.log(result[1:,0] / result[:-1,0])
    print (order)
    temp = order[order < 0.0][0]
    tempIdx = np.where(order == temp)[0][0]
    mask = (order <= 0.5)
    mask[tempIdx:] = False
    print (len(order[mask]))

    if (len(order[mask]) >= 4):
        print(result)
        raise RuntimeError("Gradient is not discrete-exact!")
    return

if __name__ == '__main__':
    config = InputParser('./optim.yml')
    optim = Optimizer(config)

    command = optim.base.forwardRunCommand('x0',True)
    fID = open('test-forward.sh','w')
    fID.write(command)
    fID.close()
    subprocess.check_call('bash test-forward.sh', shell=True)

    command = optim.base.adjointRunCommand()
    command += optim.base.gatherControlForcingGradientCommand()
    command += optim.base.innerProductCommand(optim.fl.globalGradFiles,optim.fl.globalGradFiles,optim.fl.ggFiles)
    fID = open('test-adjoint.sh','w')
    fID.write(command)
    fID.close()
    subprocess.check_call('bash test-adjoint.sh', shell=True)

    QoI0, subJ0 = optim.base.QoI('x0')
    print ("QoI0: ", QoI0)
    print ("subJ: ", subJ0)

    Grad0, subGrad0 = optim.base.readInnerProduct(optim.fl.ggFiles)
    print ("Grad0: ", Grad0)
    print ("subGrad: ", subGrad0)

    Nk = 15
    Ak = 10.0**(-4.0-0.25*np.array(range(Nk)))
    QoIk = np.zeros((Nk,),dtype=np.double)
    Gradk = np.zeros((Nk,),dtype=np.double)
    ek = np.zeros((Nk,),dtype=np.double)

    temp = ['x0/'+file for file in optim.fl.globalControlSpaceFiles]
    for k in range(optim.const.NcontrolRegion):
        temp[k] = ''

    stepFiles = optim.fl.globalControlSpaceFiles.copy()
    for k, file in enumerate(stepFiles):
        stepFiles[k] = 'x/' + file

    for k in range(Nk):
        amp = Ak[k]
        command = optim.base.zaxpyCommand(stepFiles,amp,optim.fl.globalGradFiles,temp)
        command += optim.base.forwardRunCommand('x')
        command += optim.scriptor.purgeDirectoryCommand('x')
        fID = open('test-step.sh','w')
        fID.write(command)
        fID.close()
        subprocess.check_call('bash test-step.sh', shell=True)
        QoIk[k], subJ = optim.base.QoI('x')
        Gradk[k] = (QoIk[k]-QoI0)/Ak[k]
        ek[k] = abs( (Gradk[k]-Grad0)/Grad0 )
        print ("{:.16E}".format(Ak[k]), "{:.16E}".format(QoIk[k]), "{:.16E}".format(Gradk[k]), "{:.16E}".format(ek[k]))

        fId = open(optim.fl.globalPrefix+'.gradient_accuracy.txt','a+')
        fId.write('%.16E\t%.16E\t%.16E\t%.16E\n' % (Ak[k], QoIk[k], Gradk[k], ek[k]))
        fId.close()

    checkFiniteDifferenceOrder(optim.fl.globalPrefix+'.gradient_accuracy.txt')
