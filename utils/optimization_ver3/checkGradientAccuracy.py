import numpy as np
import subprocess
from parser import parser
from inputs import InputParser
from optimizer import Optimizer, Result, Stage
import os.path as path
import h5py
#
# parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
# parser.add_argument('input_file', metavar='string', type=str,
#                     help='filename for the optimization input file.\n')

class GradientTester(Optimizer):
    fdStep = -1
    Ak = 0.0
    def __init__(self, config):
        self.fdStep = 0
        Optimizer.__init__(self, config)

        self.baseFiles = ['x0/'+file for file in self.fl.globalControlSpaceFiles]
        for k in range(self.const.NcontrolRegion):
            self.baseFiles[k] = ''

        self.stepFiles = self.fl.globalControlSpaceFiles.copy()
        for k, file in enumerate(self.stepFiles):
            self.stepFiles[k] = 'x/' + file
        return

    def printState(self):
        output = "=" * 20 + "  Gradient Tester Status  " + "=" * 20
        length = len(output)
        output += "\n"

        output += "Log file: %s\n\n" % (self.logFile)
        output += "Stage: %s\n" % (self.stage.name)
        output += "Result: %s\n" % (self.result.name)
        output += "\n"
        output += "FD step: %d\n" % (self.fdStep)
        output += "No baseline control: %s\n" % (self.zeroControlForcing)
        if (self.const.useLagrangian):
            output += "Zero Lagrangian: %s\n" % (self.zeroLagrangian)

        output += "=" * length + "\n"
        self.printAndLog(output)
        return

    def saveState(self):
        self.printState()

        with h5py.File(self.logFile, 'a') as f:
            f.attrs['fd_step'] = self.fdStep
            f.attrs['Ak'] = self.Ak
            f.attrs['control'] = self.zeroControlForcing
            f.attrs['stage'] = self.stage.value
            f.attrs['result'] = self.result.value
            if (self.const.useLagrangian):
                f.attrs['lagrangian'] = self.zeroLagrangian
        return

    def loadState(self):
        assert(path.exists(self.logFile))

        with h5py.File(self.logFile, 'r') as f:
            self.fdStep = f.attrs['fd_step']
            self.Ak = f.attrs['Ak']
            self.zeroControlForcing = f.attrs['control']
            self.stage = Stage(f.attrs['stage'])
            self.result = Result(f.attrs['result'])
            if (self.const.useLagrangian):
                self.zeroLagrangian = f.attrs['lagrangian']
        return

    def schedule(self):
        if ((self.stage is Stage.CG_FW) and (self.result is Result.UNEXECUTED)):
            command = self.base.forwardRunCommand('x0',True)
            self.writeCommandFile(command, self.fl.globalCommandFile)
            self.saveState()
            return

        assert(self.result is Result.SUCCESS)
        # Based on the completed stage, determine next action
        if (self.stage is Stage.CG_FW):
            command = self.base.adjointRunCommand()
            command += self.base.gatherControlForcingGradientCommand()
            command += self.base.innerProductCommand(self.fl.globalGradFiles,
                                    self.fl.globalGradFiles, self.fl.ggFiles)
            self.writeCommandFile(command, self.fl.globalCommandFile)
            self.stage, self.result = Stage.CG_AD, Result.UNEXECUTED
        elif (self.stage is Stage.CG_AD):
            QoI0, subJ0 = self.base.QoI('x0')
            self.printAndLog("QoI0: %.15E" % QoI0)
            fwDset = "gradient_test/objectives"
            self.appendDataset(fwDset, np.array([[QoI0]+list(subJ0)]))

            Grad0, subGrad0 = self.base.readInnerProduct(optim.fl.ggFiles)
            self.printAndLog("Grad0: %.15E" % Grad0)
            adjDset = "gradient_test/gradients"
            self.appendDataset(adjDset, np.array([[Grad0]+list(subGrad0)]))

            self.Ak = 10.0 ** (-4.0 - 0.25 * self.fdStep)
            command = self.base.zaxpyCommand(self.stepFiles,self.Ak,self.fl.globalGradFiles,self.baseFiles)
            command += self.base.forwardRunCommand('x')
            command += self.scriptor.purgeDirectoryCommand('x')
            self.writeCommandFile(command, self.fl.globalCommandFile)

            self.stage, self.result = Stage.BRKT, Result.UNEXECUTED
        elif (self.stage is Stage.BRKT):
            with h5py.File(self.logFile, 'r') as f:
                QoI0 = f['gradient_test/objectives'][0,0]
                Grad0 = f['gradient_test/gradients'][0,0]
            QoIk, subJ = optim.base.QoI('x')
            Gradk = (QoIk-QoI0) / self.Ak
            ek = abs( (Gradk-Grad0)/Grad0 )
            self.printAndLog("%.15E\t%.15E\t%.15E\t%.15E" % (self.Ak, QoIk, Gradk, ek))
            self.appendDataset('gradient_test/accuracy', np.array([[self.Ak, QoIk, Gradk, ek]]))

            self.fdStep += 1
            if (self.fdStep == 15):
                self.checkFiniteDifferenceOrder()
                self.writeCommandFile('exit 1\n', self.fl.globalCommandFile)
                return

            self.Ak = 10.0 ** (-4.0 - 0.25 * self.fdStep)
            command = self.base.zaxpyCommand(self.stepFiles,self.Ak,self.fl.globalGradFiles,self.baseFiles)
            command += self.base.forwardRunCommand('x')
            command += self.scriptor.purgeDirectoryCommand('x')
            self.writeCommandFile(command, self.fl.globalCommandFile)

            self.stage, self.result = Stage.BRKT, Result.UNEXECUTED

        self.saveState()
        return

    def checkFiniteDifferenceOrder(self):
        with h5py.File(self.logFile, 'r') as f:
            result = f['gradient_test/accuracy'][...]
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
    args = parser.parse_args()
    if ((args.mode == "log_result") and (args.result is None)):
        raise RuntimeError("Specify the result as '--result integer'!\n")

    config = InputParser(args.input_file)
    optim = GradientTester(config)

    if (args.mode == "schedule"):
        optim.schedule()
    elif (args.mode == "log_result"):
        optim.checkResult(args.result)
    elif (args.mode == "setup"):
        optim.setupDirectories()

    else:
        print("choose the action. for possible choices, add '-h' flag and run it again.")
    # args = parser.parse_args()
    # config = InputParser(args.input_file)
    # optim = Optimizer(config)
    #
    # command = optim.base.forwardRunCommand('x0',True)
    # fID = open('test-forward.sh','w')
    # fID.write(command)
    # fID.close()
    # subprocess.check_call('bash test-forward.sh', shell=True)
    #
    # command = optim.base.adjointRunCommand()
    # command += optim.base.gatherControlForcingGradientCommand()
    # command += optim.base.innerProductCommand(optim.fl.globalGradFiles,optim.fl.globalGradFiles,optim.fl.ggFiles)
    # fID = open('test-adjoint.sh','w')
    # fID.write(command)
    # fID.close()
    # subprocess.check_call('bash test-adjoint.sh', shell=True)
    #
    # QoI0, subJ0 = optim.base.QoI('x0')
    # print ("QoI0: ", QoI0)
    # print ("subJ: ", subJ0)
    #
    # Grad0, subGrad0 = optim.base.readInnerProduct(optim.fl.ggFiles)
    # print ("Grad0: ", Grad0)
    # print ("subGrad: ", subGrad0)
    #
    # Nk = 15
    # Ak = 10.0**(-4.0-0.25*np.array(range(Nk)))
    # QoIk = np.zeros((Nk,),dtype=np.double)
    # Gradk = np.zeros((Nk,),dtype=np.double)
    # ek = np.zeros((Nk,),dtype=np.double)
    #
    # temp = ['x0/'+file for file in optim.fl.globalControlSpaceFiles]
    # for k in range(optim.const.NcontrolRegion):
    #     temp[k] = ''
    #
    # stepFiles = optim.fl.globalControlSpaceFiles.copy()
    # for k, file in enumerate(stepFiles):
    #     stepFiles[k] = 'x/' + file
    #
    # for k in range(Nk):
    #     amp = Ak[k]
    #     command = optim.base.zaxpyCommand(stepFiles,amp,optim.fl.globalGradFiles,temp)
    #     command += optim.base.forwardRunCommand('x')
    #     command += optim.scriptor.purgeDirectoryCommand('x')
    #     fID = open('test-step.sh','w')
    #     fID.write(command)
    #     fID.close()
    #     subprocess.check_call('bash test-step.sh', shell=True)
    #     QoIk[k], subJ = optim.base.QoI('x')
    #     Gradk[k] = (QoIk[k]-QoI0)/Ak[k]
    #     ek[k] = abs( (Gradk[k]-Grad0)/Grad0 )
    #     print ("{:.16E}".format(Ak[k]), "{:.16E}".format(QoIk[k]), "{:.16E}".format(Gradk[k]), "{:.16E}".format(ek[k]))
    #
    #     fId = open(optim.fl.globalPrefix+'.gradient_accuracy.txt','a+')
    #     fId.write('%.16E\t%.16E\t%.16E\t%.16E\n' % (Ak[k], QoIk[k], Gradk[k], ek[k]))
    #     fId.close()
    #
    # checkFiniteDifferenceOrder(optim.fl.globalPrefix+'.gradient_accuracy.txt')
