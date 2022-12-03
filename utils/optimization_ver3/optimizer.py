from command_scriptor import scriptorSwitcher
from constants import Constants
from filenames import FilenameList
from base import BaseCommander
from base_extension import BaseCommanderExtended
from enum import Enum
import numpy as np
import subprocess
import pandas as pd
import os
import os.path as path
import h5py
import logging

__all__ = ['Optimizer']

magudiExes = ['forward', 'adjoint', 'control_space_norm', 'qfile_zaxpy',
              'spatial_inner_product', 'zxdoty', 'zwxmwy', 'zaxpy',
              'slice_control_forcing', 'paste_control_forcing', 'patchup_qfile']

class Result(Enum):
    UNEXECUTED = -1
    SUCCESS = 0
    FAIL = 1

class Stage(Enum):
    CG_FW = 0
    CG_AD = 1
    POST_AD = 2
    BRKT = 3
    LNMN = 4

class Optimizer:
    config = None
    scriptor = None
    const = None
    fl = None
    base = None

    binaryPath = ''

    inputFile = ''
    logFile = ''
    journalFile = ''

    isInitial = True   # the very beginning of the entire multi-point optimization
    hyperStep = -1     # Penalty-method level hyper-iteration
    cgStep = -1        # Number of CG step at the current hyper step
    lineStep = -1      # Number of line step at the current CG step
    bracketed = False
    linminInitial = True

    stage = None                # Scheduled action at the current (hyperStep, cgStep, lineStep)
    result = Result.UNEXECUTED  # Result of scheduled action at the current (hyperStep, cgStep, lineStep)

    zeroControlForcing = True   # control at the current CG step (directory x0)
    zeroLagrangian = True       # augmented_lagrangian at the current stage (root directory)

    def __init__(self, config):
        self.config = config
        self.const = Constants(config)
        self.fl = FilenameList(config, self.const)

        scriptorType = config.getInput(['command_scriptor', 'type'], datatype=str)
        self.scriptor = scriptorSwitcher[scriptorType](config, self.fl)

        if (config.getInput(['time_splitting', 'use_state_mollifier'], fallback=False)):
            self.base = BaseCommanderExtended(config, self.scriptor, self.const, self.fl)
        else:
            self.base = BaseCommander(config, self.scriptor, self.const, self.fl)

        self.binaryPath = self.config.getInput(['magudi', 'binary_directory'], datatype=str)
        self.binaryPath = os.path.abspath(self.binaryPath)

        defaultJournal = "%s.journal.txt" % self.fl.globalPrefix
        self.journalFile = config.getInput(['optimization', 'journal_file'], fallback=defaultJournal)
        logging.basicConfig(filename=self.journalFile, filemode='a',
                            level=logging.INFO, format='%(asctime)s\n%(message)s')

        self.inputFile = config.filename
        self.logFile = config.getInput(['optimization', 'log_file'], datatype=str)
        assert(self.logFile[-3:] == ".h5")
        self.isInitial = (not path.exists(self.logFile))

        if (self.isInitial):
            self.hyperStep = 0
            self.cgStep = 0
            self.lineStep = 0
            self.stage, self.result = Stage.CG_FW, Result.UNEXECUTED
            self.bracketed = False
            self.linminInitial = True
            self.zeroControlForcing = config.getInput(['optimization', 'initial', 'zero_control'], fallback=True)
            if (self.const.useLagrangian):
                self.zeroLagrangian = config.getInput(['optimization', 'initial', 'zero_lagrangian'], fallback=True)

            self.saveState()

            with h5py.File(self.logFile, 'a') as f:
                dset = f.create_dataset("%s/columns" % self.fl.forwardLog,
                                        data=list(self.fl.forwardLogColumns.keys()))
                dset = f.create_dataset("%s/columns" % self.fl.gradientLog,
                                        data=list(self.fl.gradientLogColumns.keys()))
                dset = f.create_dataset("%s/columns" % self.fl.lineMinLog,
                                        data=self.fl.lineMinLogColumns)
        else:
            self.loadState()
        return

    def schedule(self):
        if ((self.stage is Stage.CG_FW) and (self.result is Result.UNEXECUTED)):
            self.initialForward()
            self.saveState()
            return

        assert(self.result is Result.SUCCESS)
        # Based on the completed stage, determine next action
        if (self.stage is Stage.CG_FW):
            if (self.cgStep == 0):
                self.initialAdjoint()
            else:
                case = self.afterLinmin()
            self.stage, self.result = Stage.CG_AD, Result.UNEXECUTED
        elif (self.stage is Stage.CG_AD):
            success = self.beforeLinmin()
            if (success!=0):
                self.printAndLog('pre-processing for line minimization is failed.')
            self.stage, self.result = Stage.POST_AD, Result.UNEXECUTED
        elif (self.stage is Stage.POST_AD):
            success = self.setupInitialSteps()
            if (success != 0):
                self.printAndLog('initialization for mnbrak is failed.')
            self.lineStep += 1
            self.stage, self.result = Stage.BRKT, Result.UNEXECUTED
        elif (self.stage is Stage.BRKT):
            self.bracketed = self.nextMnbrak()
            if (self.bracketed):
                self.stage, self.linminInitial = Stage.LNMN, True
                self.schedule()
                return
            else:
                self.lineStep += 1
                self.stage, self.result = Stage.BRKT, Result.UNEXECUTED
        elif (self.stage is Stage.LNMN):
            continueLinmin = self.nextLinmin(self.linminInitial)#, stop=linminStop)
            self.linminInitial = False
            if (continueLinmin):
                self.lineStep += 1
                self.stage, self.result = Stage.LNMN, Result.UNEXECUTED
            else:
                self.cgStep += 1
                self.lineStep = 0
                self.zeroControlForcing, self.bracketed = False, False
                self.stage, self.result = Stage.CG_FW, Result.UNEXECUTED

        self.saveState()
        return

    def checkResult(self, result):
        self.result = Result(result) if (result == 0) else Result.FAIL
        self.saveState()

        if (self.result is Result.FAIL):
            raise RuntimeError("%s is not run successfully." % self.fl.globalCommandFile)
        return

    def loadState(self):
        assert(path.exists(self.logFile))

        with h5py.File(self.logFile, 'r') as f:
            self.hyperStep = f.attrs['hyper_step']
            self.cgStep = f.attrs['cg_step']
            self.lineStep = f.attrs['line_step']
            self.bracketed = f.attrs['bracketed']
            self.linminInitial = f.attrs['linmin_initial']
            self.zeroControlForcing = f.attrs['control']
            self.stage = Stage(f.attrs['stage'])
            self.result = Result(f.attrs['result'])
            if (self.const.useLagrangian):
                self.zeroLagrangian = f.attrs['lagrangian']
        return

    def loadPreviousPenalty(self):
        assert(path.exists(self.logFile))
        assert(self.const.useLagrangian)
        assert((self.hyperStep > 0) and (self.cgStep == 0) and (self.lineStep == 0))

        with h5py.File(self.logFile, 'r') as f:
            dsetName = "%d/penalty_weight" % (self.hyperStep - 1)
            weight = f[dsetName][...]
        return weight

    def appendDataset(self, dataset, data):
        with h5py.File(self.logFile, 'a') as f:
            if (dataset in f):
                self.printAndLog("Appending %s.." % dataset)
                f[dataset].resize(f[dataset].shape[0] + 1, axis=0)
                f[dataset][-1] = data
            else:
                self.printAndLog("Creating %s.." % dataset)
                dset = f.create_dataset(dataset, data=data, maxshape=(None, data.shape[1]))
        self.printAndLog("Succeeded.")
        return

    def saveCGStepData(self):
        self.printAndLog("Saving data at cg step %d/%d.." % (self.hyperStep, self.cgStep))
        fwHist = np.empty(shape=(0,4))
        adjHist = np.empty(shape=(0,4))
        for k in range(self.const.Nsplit):
            fwHist = np.append(fwHist, np.loadtxt("x0/%s" % (self.fl.forwardHistoryFiles[k])), axis=0)
            adjHist = np.append(adjHist, np.loadtxt("x0/%s" % self.fl.gradientHistoryFiles[k]), axis=0)

        if (self.cgStep > 0):
            with h5py.File(self.logFile, 'a') as f:
                dsetName = "%s/%d/%d" % (self.fl.forwardHistory, self.hyperStep, self.cgStep)
                dset = f.create_dataset(dsetName, data=fwHist)
                dsetName = "%s/%d/%d" % (self.fl.gradientHistory, self.hyperStep, self.cgStep)
                dset = f.create_dataset(dsetName, data=adjHist)

                dsetName = "%s/%d/%d" % (self.fl.lineMinLog, self.hyperStep, self.cgStep)
                df = pd.read_csv(self.fl.lineMinLogFile, sep='\t', header=0)
                data = np.array(df.iloc[:, :2])
                dset = f.create_dataset(dsetName, data=data)
                for dir in ['a', 'b', 'c']:
                    idx = np.array(df['step'][df['directory index']==dir])[0]
                    dset.attrs[dir] = idx
        self.printAndLog("Succeeded.")
        return

    def saveState(self):
        self.printState()

        with h5py.File(self.logFile, 'a') as f:
            f.attrs['hyper_step'] = self.hyperStep
            f.attrs['cg_step'] = self.cgStep
            f.attrs['line_step'] = self.lineStep
            f.attrs['bracketed'] = self.bracketed
            f.attrs['linmin_initial'] = self.linminInitial
            f.attrs['control'] = self.zeroControlForcing
            f.attrs['stage'] = self.stage.value
            f.attrs['result'] = self.result.value
            if (self.const.useLagrangian):
                f.attrs['lagrangian'] = self.zeroLagrangian
        return

    def printState(self):
        output = "=" * 20 + "  Optimizer Status  " + "=" * 20
        length = len(output)
        output += "\n"

        output += "Log file: %s\n\n" % (self.logFile)
        output += "Stage: %s\n" % (self.stage.name)
        output += "Result: %s\n" % (self.result.name)
        output += "\n"
        output += "Hyper step: %d\n" % (self.hyperStep)
        output += "CG step: %d\n" % (self.cgStep)
        output += "Line step: %d\n" % (self.lineStep)
        output += "Bracketed: %s\n" % (self.bracketed)
        output += "LinMin Initial: %s\n" % (self.linminInitial)
        output += "No baseline control: %s\n" % (self.zeroControlForcing)
        if (self.const.useLagrangian):
            output += "Zero Lagrangian: %s\n" % (self.zeroLagrangian)

        output += "=" * length + "\n"
        self.printAndLog(output)
        return

    def printAndLog(self, message):
        print(message)
        logging.info(message)
        return

    def writeCommandFile(self, command, filename):
        fID = open(filename, 'w')
        fID.write(command)
        fID.close()
        return

    def parabolic_interp(self, stepBracket, JBracket):
        a,b,c = stepBracket
        fa,fb,fc = JBracket

        x = b - 0.5 * ( ((b-a)**2)*(fb-fc) - ((b-c)**2)*(fb-fa) )                                           \
                    / ( (b-a)*(fb-fc) - (b-c)*(fb-fa) )
        return x

    def setMagudiInputFiles(self, inputFiles, key, vals):
        assert(len(inputFiles) == len(vals))
        commandString = self.base.generalSetOptionCommand

        commands = []
        for inputFile, val in zip(inputFiles, vals):
            commands += ['./setOption.sh %s "%s" %s ' % (inputFile, key, val)]
        commandString += self.scriptor.nonMPILoopCommand(commands, key)
        return commandString

    def setupDirectories(self):
        prerequisites = self.fl.globalNormFiles
        for k in range(self.const.Nsplit):
            idx = self.const.startTimestep + k * self.const.Nts
            prerequisites += ["%s-%08d.q" % (self.fl.globalPrefix, idx)]
        for file in prerequisites:
            assert(path.exists(file))

        commandString = ''

        for exe in magudiExes:
            commandString += 'ln -s %s/%s ./\n' % (self.binaryPath, exe)

        for dir in self.fl.dirList:
            if (not path.exists(dir)):
                commandString += 'mkdir -p %s\n' % dir

        commandString +=                                                                \
        'for dir in {a,b,c,x,x0}\n'                                                     \
        'do\n'                                                                          \
        '   mkdir $dir\n'                                                               \
        '   for k in {0..%d}\n' % (self.const.Nsplit - 1)
        commandString +=                                                                \
        '   do\n'                                                                       \
        '       mkdir -p ${dir}/${k}\n'                                                 \
        '       ln -s %s/forward ./${dir}/${k}/\n' % self.binaryPath
        commandString +=                                                                \
        '       ln -s %s/adjoint ./${dir}/${k}/\n' % self.binaryPath
        commandString +=                                                                \
        '       cp ./%s ./${dir}/${k}/magudi-${k}.inp\n' % self.fl.globalInputFile
        commandString +=                                                                \
        '   done\n'                                                                     \
        'done\n'
        commandString += 'if [ $? -ne 0 ]; then exit -1; fi\n\n'

        targetInputFiles = []
        for bdir in ['a', 'b', 'c', 'x', 'x0']:
            targetInputFiles += ['%s/%s/%s'%(bdir,dir,file) for dir, file in zip(self.fl.directories, self.fl.inputFiles)]

        rootFiles = self.config.getInput(['magudi', 'root_files'], datatype=dict)
        for key, val in rootFiles.items():
            assert(type(val) is str)
            vals = ['\'"../../%s"\'' % val] * len(targetInputFiles)
            commandString += self.setMagudiInputFiles(targetInputFiles, key, vals)

        vals = ['\'"%s-%d"\'' % (self.fl.globalPrefix, k) for k in range(self.const.Nsplit)] * 5
        commandString += self.setMagudiInputFiles(targetInputFiles, 'output_prefix', vals)

        prefix = self.fl.globalPrefix
        if (type(self.base) is not BaseCommanderExtended): prefix = '../' + prefix
        vals = ['\'"%s-%d.ic.q"\'' % (prefix, k) for k in range(self.const.Nsplit)] * 5
        commandString += self.setMagudiInputFiles(targetInputFiles, 'initial_condition_file', vals)

        vals = ['%d' % self.const.Nts] * len(targetInputFiles)
        commandString += self.setMagudiInputFiles(targetInputFiles, 'number_of_timesteps', vals)

        vals = ['%d' % self.const.Nts] * (4 * self.const.Nsplit)
        adjointSaveTimestep = self.config.getInput(['magudi', 'adjoint_save_timesteps'], fallback = self.const.Nts)
        vals += ['%d' % adjointSaveTimestep] * self.const.Nsplit
        commandString += self.setMagudiInputFiles(targetInputFiles, 'save_interval', vals)

        vals = ['"true"'] * (4 * self.const.Nsplit)
        vals += ['"%s"' % ('true' if (not self.zeroControlForcing) else 'false')] * self.const.Nsplit
        commandString += self.setMagudiInputFiles(targetInputFiles, 'controller_switch', vals)

        commands = []
        for k in range(self.const.Nsplit):
            idx = self.const.startTimestep + k * self.const.Nts
            oldFile = "%s-%08d.q" % (self.fl.globalPrefix, idx)
            newFile = "%s-%d.ic.q" % (self.fl.globalPrefix, k)
            commands += ['cp %s x0/%s' % (oldFile, newFile)]
        commandString += self.scriptor.nonMPILoopCommand(commands, 'copying_initial_conditions')

        self.writeCommandFile(commandString, self.fl.globalCommandFile)
        return

    def initialForward(self):
        command = ''
        command += self.base.forwardRunCommand('x0', self.zeroControlForcing)

        if (self.const.useLagrangian and (self.hyperStep > 0)):
            initialLagrangian = (not self.zeroLagrangian)
            weight = self.loadPreviousPenalty()
            command += optim.base.updateLagrangian(weight, initialLagrangian)
        self.writeCommandFile(command, self.fl.globalCommandFile)
        return

    def initialAdjoint(self):
        command = ''
        command += self.base.adjointRunCommand()
        command += self.base.gatherControlForcingGradientCommand()
        command += self.base.innerProductCommand(self.fl.globalGradFiles,
                                                 self.fl.globalGradFiles,
                                                 self.fl.ggFiles)
        self.writeCommandFile(command, self.fl.globalCommandFile)
        return

    def setupInitialSteps(self):
        if(self.const.saveStateLog):
            self.readStateDistance()

        commandString = ''
        J0, dummy = self.base.QoI()

        commands = []
        for k in range(1, self.const.Nsplit):
            commands += ['cp x0/%s ./a/ ' % self.fl.icFiles[k]]
        if (not self.zeroControlForcing):
            for j in range(self.const.NcontrolRegion):
                commands += ['cp x0/%s ./a/ ' % self.const.globalControlSpaceFiles[j]]
        commandString += self.scriptor.nonMPILoopCommand(commands,'copy_control_params')

        numFiles = len(os.listdir(self.fl.LINMINDIR))
        if (numFiles>0):
            df_file = '%s/%d/%s' % (self.fl.LINMINDIR, numFiles-1, self.fl.lineMinLogFile)
            df_temp = pd.read_csv(df_file, sep='\t', header=0)
            lastStep = np.array(df_temp['step'][df_temp['directory index']=='b'])[0]

        steps, Js = np.zeros(2), np.zeros(2)
        steps[1] = self.const.initial_step if numFiles==0 else lastStep
        Js[0] = J0

        temp = ['x0/'+file for file in self.fl.globalControlSpaceFiles]
        if (self.zeroControlForcing):
            for j in range(self.const.NcontrolRegion):
                temp[j] = ''
        target = self.fl.globalControlSpaceFiles.copy()
        for k, file in enumerate(target):
            target[k] = 'x/' + file
        commandString += self.base.zaxpyCommand(target, -steps[1], self.fl.globalConjugateGradientFiles, temp)
        commandString += '\n'
        commandString += self.base.forwardRunCommand('x')
        commandString += self.scriptor.purgeDirectoryCommand('x')
        self.writeCommandFile(commandString, self.fl.globalCommandFile)

        data = {'step':steps,'QoI':Js,'directory index':['a','x']}
        df = pd.DataFrame(data)
        df.to_csv(self.fl.lineMinLogFile, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)

        self.printAndLog('Initial steps written in command. Run ' + self.fl.globalCommandFile + '.')
        return 0

    def nextMnbrak(self):
        df = self.base.collectQoI(self.fl.lineMinLogFile)
        steps = np.array(df['step'][df['directory index']!='0'])
        Js = np.array(df['step'][df['directory index']!='0'])
        dirIdx = np.array(df['directory index'][df['directory index']!='0'])
        N = steps.size

        a = np.array(df['step'][df['directory index']=='a'])[0]
        Ja = np.array(df['QoI'][df['directory index']=='a'])[0]
        x = np.array(df['step'][df['directory index']=='x'])[0]
        Jx = np.array(df['QoI'][df['directory index']=='x'])[0]
        bExists = sum(dirIdx=='b')
        cExists = sum(dirIdx=='c')

        if (Jx > self.const.huge):
            self.printAndLog('Forward run crashed. Retrying with a smaller step.')
            df.loc[df['directory index']=='x','directory index'] = '0'
            Cr = 1.0 - 1.0 / self.const.golden_ratio
            if (bExists):
                b = np.array(df['step'][df['directory index']=='b'])[0]
                new_x = b + Cr * (x-b)
            else:
                new_x = a + Cr * (x-a)
        elif (bExists):
            b = np.array(df['step'][df['directory index']=='b'])[0]
            Jb = np.array(df['QoI'][df['directory index']=='b'])[0]

            if (Jx > Jb):
                self.base.switchDirectory('x','c',df)

                df.to_csv(self.fl.lineMinLogFile, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)
                self.printAndLog(np.array(df[df['directory index']!='0']))
                self.printAndLog('MNBRAK: initial mininum bracket is prepared.')
                return True # bracketed.
            else:
                self.base.switchDirectory('x','b',df)
                self.base.switchDirectory('a','x',df)
                df.loc[df['directory index']=='x','directory index'] = '0'

                if (len(df[(df['QoI'] > self.const.huge)]) > 0):
                    maxX = np.min(np.array(df['step'][df['QoI'] > self.const.huge]))
                    Cr = 1.0 - 1.0 / self.const.golden_ratio
                    new_x = x + Cr * (maxX - x)
                else:
                    new_x = x * self.const.golden_ratio if (x < self.const.safe_zone) else x * 1.05
                self.printAndLog('MNBRAK: expanding the bracket - Run intermediate forward simulation.')
        elif (cExists):
            c = np.array(df['step'][df['directory index']=='c'])[0]
            Jc = np.array(df['QoI'][df['directory index']=='c'])[0]

            if (Jx < Ja):
                self.base.switchDirectory('x','b',df)

                df.to_csv(self.fl.lineMinLogFile, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)
                self.printAndLog(np.array(df[df['directory index']!='0']))
                self.printAndLog('MNBRAK: initial mininum bracket is prepared.')
                return True # bracketed.
            else:
                self.base.switchDirectory('x','c',df)
                df.loc[df['directory index']=='x','directory index'] = '0'

                new_x = x / self.const.golden_ratio
                self.printAndLog('MNBRAK: narrowing the bracket - Run intermediate forward simulations.')
        else:
            if (Jx > Ja):
                self.base.switchDirectory('x','c',df)

                new_x = x / self.const.golden_ratio
                self.printAndLog('MNBRAK: narrowing the bracket - Run intermediate forward simulations.')
            else:
                self.base.switchDirectory('x','b',df)

                if (len(df[(df['QoI'] > self.const.huge)]) > 0):
                    maxX = np.min(np.array(df['step'][df['QoI'] > self.const.huge]))
                    Cr = 1.0 - 1.0 / self.const.golden_ratio
                    new_x = x + Cr * (maxX-x)
                else:
                    new_x = x * self.const.golden_ratio if (x < self.const.safe_zone) else x * 1.05
                self.printAndLog('MNBRAK: expanding the bracket - Run intermediate forward simulation.')

        data = {'step':[new_x], 'QoI':[np.nan], 'directory index':['x']}
        new_df = pd.DataFrame(data)
        df = df.append(new_df, ignore_index=True)

        commandString = ''
        temp = ['x0/'+file for file in self.fl.globalControlSpaceFiles]
        if (self.zeroControlForcing):
            for j in range(self.const.NcontrolRegion):
                temp[j] = ''
        target = ['x/'+file for file in self.fl.globalControlSpaceFiles]
        commandString += self.scriptor.parallelPurgeCommand(target,'purge_target')
        commandString += self.base.zaxpyCommand(target, -new_x, self.fl.globalConjugateGradientFiles, temp)
        commandString += '\n'
        commandString += self.base.forwardRunCommand('x')
        commandString += self.scriptor.purgeDirectoryCommand('x')
        self.writeCommandFile(commandString, self.fl.globalCommandFile)

        df.to_csv(self.fl.lineMinLogFile, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)
        self.printAndLog(df[df['directory index']!='0'])

        return False # not bracketed.

    def nextLinmin(self, initial=True, stop=False):
        if (initial):
            df = pd.read_csv(self.fl.lineMinLogFile, sep='\t', header=0)
        else:
            df = self.base.collectQoI(self.fl.lineMinLogFile)

            steps = np.array(df['step'][df['directory index'] != '0'])
            Js = np.array(df['QoI'][df['directory index'] != '0'])
            dirIdx = np.array(df['directory index'][df['directory index'] != '0'])
            N = steps.size
            if (N != 4):
                raise RuntimeError('Cannot find values corresponding to a, b, c, and x!')
            minIdx = max(np.argmin(Js), 1)

            temp = ['a','b','c']
            for k in range(3):
                self.base.switchDirectory(dirIdx[minIdx-1+k], temp[k], df)
                dirIdx = np.array(df['directory index'][df['directory index'] != '0'])
            df.loc[df['directory index']=='x','directory index'] = '0'
            df.to_csv(self.fl.lineMinLogFile, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)

        if (stop):
            self.printAndLog(df[df['directory index']!='0'])
            self.printAndLog('LINMIN: line minimization is stopped.')

            self.writeCommandFile('exit 0\n', self.fl.globalCommandFile)
            return 0

        a = np.array(df['step'][df['directory index']=='a'])[0]
        b = np.array(df['step'][df['directory index']=='b'])[0]
        c = np.array(df['step'][df['directory index']=='c'])[0]
        Ja = np.array(df['QoI'][df['directory index']=='a'])[0]
        Jb = np.array(df['QoI'][df['directory index']=='b'])[0]
        Jc = np.array(df['QoI'][df['directory index']=='c'])[0]

        if ( (c-a) < b * self.const.linminTol + self.const.eps ):
            self.printAndLog('steps: %.16E %.16E %.16E' % (a,b,c))
            self.printAndLog('QoIs: %.16E %.16E %.16E' % (Ja,Jb,Jc))
            self.printAndLog('LINMIN: line minimization is finished.')

            commandString = ''

            target = ['a/'+file for file in self.fl.globalControlSpaceFiles]
            commandString += self.scriptor.parallelPurgeCommand(target,'purge_a')
            target = ['c/'+file for file in self.fl.globalControlSpaceFiles]
            commandString += self.scriptor.parallelPurgeCommand(target,'purge_c')
            target = ['x/'+file for file in self.fl.globalControlSpaceFiles]
            commandString += self.scriptor.parallelPurgeCommand(target,'purge_x')

            commands = []
            for k in range(self.const.NcontrolSpace):
                commands += ['cp b/%s x0/ ' % (self.fl.globalControlSpaceFiles[k])]
                commands += ['mv %s %s' % (self.fl.globalGradFiles[k], self.fl.previousGradFiles[k])]
                commands += ['mv %s %s' % (self.fl.globalConjugateGradientFiles[k], self.fl.previousCGFiles[k])]
            commandString += self.scriptor.nonMPILoopCommand(commands,'move_line_minimum_files')

            commandString += '\n'
            if( self.zeroControlForcing and (not self.const.ignoreController) ):
                commandString += self.base.generalSetOptionCommand
                targetInputFiles = ['x0/%s/%s'%(dir,file) for dir, file in zip(self.fl.directories, self.fl.inputFiles)]
                commandString += self.setMagudiInputFiles(targetInputFiles, "controller_switch", ["true"] * self.const.Nsplit)

            commandString += self.base.forwardRunCommand()

            self.writeCommandFile(commandString, self.fl.globalCommandFile)
            return False # completed linmin.

        xs = self.parabolic_interp([a,b,c],[Ja,Jb,Jc])
        Cr = 1.0 - 1.0 / self.const.golden_ratio
        if( (xs>c) or (xs<a) or (abs(np.log10((c-b)/(b-a)))>1.0) ):
            if( b>0.5*(a+c) ):
                xs = b - Cr * (b-a)
            else:
                xs = b + Cr * (c-b)
        data = {'step':[xs], 'QoI':[np.nan], 'directory index':['x']}
        new_df = pd.DataFrame(data)
        df = df.append(new_df, ignore_index=True)

        commandString = ''
        temp = ['x0/'+file for file in self.fl.globalControlSpaceFiles]
        if (self.zeroControlForcing):
            for j in range(self.const.NcontrolRegion):
                temp[j] = ''
        target = ['x/'+file for file in self.fl.globalControlSpaceFiles]
        commandString += self.scriptor.parallelPurgeCommand(target,'purge_target')
        commandString += self.base.zaxpyCommand(target, -xs, self.fl.globalConjugateGradientFiles, temp)
        commandString += '\n'
        commandString += self.base.forwardRunCommand('x')
        commandString += self.scriptor.purgeDirectoryCommand('x')
        self.writeCommandFile(commandString, self.fl.globalCommandFile)

        df.to_csv(self.fl.lineMinLogFile, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)
        self.printAndLog(df[df['directory index']!='0'])
        self.printAndLog('LINMIN: next linmin evaluation is prepared-')
        return True # continue linmin.

    def beforeLinmin(self):
        J0, subJ0 = self.base.QoI()
        gg, subgg = self.base.readInnerProduct(self.fl.ggFiles)
        self.saveCGStepData()

        fwDset = "%s/%d" % (self.fl.forwardLog, self.hyperStep)
        adjDset = "%s/%d" % (self.fl.gradientLog, self.hyperStep)
        if (self.cgStep > 0):
            with h5py.File(self.logFile, 'r') as f:
                J1 = f[fwDset][-1, self.fl.forwardLogColumns['total']]
                gg1 = f[adjDset][-1, self.fl.gradientLogColumns['total']]
                gg0 = f[adjDset][0, self.fl.gradientLogColumns['total']]

        self.appendDataset(fwDset, np.array([[J0]+list(subJ0)]))
        self.appendDataset(adjDset, np.array([[gg]+list(subgg)]))

        commandString = ''
        if (self.const.saveStateLog):
            commandString += self.base.checkStateDistance()

        # Initial cg step
        if (self.cgStep == 0):
            commands = []
            for k in range(self.const.NcontrolSpace):
                commands += ['cp %s %s' % (self.fl.globalGradFiles[k], self.fl.globalConjugateGradientFiles[k])]
            commandString += self.scriptor.nonMPILoopCommand(commands,'initial-conjugate-gradient')
            self.writeCommandFile(commandString, self.fl.globalCommandFile)

            self.printAndLog('Initial line minimization is ready. Run mnbrak and linmin procedures.')
            return 0

        #reduction = gg/gg0
        #if (reduction<=tol):
        if (gg <= self.const.tol):
            self.printAndLog('FRPRMN - after linmin: conjugate-gradient optimization is finished.')
            self.writeCommandFile('exit 1\n', self.fl.globalCommandFile)
            return 0

        dgg, dummy = self.base.readInnerProduct(self.fl.dggFiles)

        # Fletcher-Reeves
        gamma1 = gg/gg1
        # Polak-Ribiere
        gamma = dgg/gg1
        if (gamma > gamma1):
            gamma = gamma1
        elif (gamma < -gamma1):
            gamma = -gamma1

        commandString += self.base.zaxpyCommand(self.fl.globalConjugateGradientFiles, gamma,
                                                self.base.previousCGFiles, self.base.globalGradFiles)
        commandString += '\n'
        commandString += self.scriptor.parallelPurgeCommand(self.fl.previousCGFiles,'purge_prev_cg')
        self.writeCommandFile(commandString, self.fl.globalCommandFile)

        self.printAndLog('line minimization is ready. Run mnbrak and linmin procedures.')
        return 0

    def afterLinmin(self):
        commandString = ''

        commandString += '\n'
        commandString += self.base.adjointRunCommand()
        commandString += '\n'
        commandString += self.base.gatherControlForcingGradientCommand()
        commandString += '\n'
        commandString += self.base.innerProductCommand(self.fl.globalGradFiles, self.fl.globalGradFiles, self.fl.ggFiles)
        commandString += self.scriptor.purgeDirectoryCommand('x0')

        if (self.const.saveDiffFiles):
            subprocess.check_call('mkdir -p diffLog/%d'%numFiles,shell=True)
            for k in range(self.const.Nsplit - 1, self.const.Nsplit - self.const.Ndiff * self.const.diffStep, -self.const.diffStep):
                commands += ['mv %s diffLog/%d/%s-%d.diff.%d.q' % (self.fl.diffFiles[k],numFiles,self.fl.globalPrefix,k,numFiles)]

        # Polak-Ribiere
        commandString += self.base.dggCommand()
        commandString += self.scriptor.parallelPurgeCommand(self.fl.previousGradFiles,'purge_prev_grad')
        self.writeCommandFile(commandString, self.fl.globalCommandFile)

        self.printAndLog('FRPRMN - after linmin: postprocessing is finished. Run new forward/adjoint simulations.')
        return 1
