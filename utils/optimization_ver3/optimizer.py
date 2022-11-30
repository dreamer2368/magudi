from command_scriptor import scriptorSwitcher
from constants import Constants
from filenames import FilenameList
from base import BaseCommander
from base_extension import BaseCommanderExtended
import numpy as np
import subprocess
import pandas as pd
import os

__all__ = ['Optimizer']

class Optimizer:
    scriptor = None
    const = None
    fl = None
    base = None

    inputFile = ''

    def __init__(self, config):
        self.const = Constants(config)
        self.fl = FilenameList(config, self.const)

        scriptorType = config.getInput(['command_scriptor', 'type'], datatype=str)
        self.scriptor = scriptorSwitcher[scriptorType](config, self.fl)

        if (config.getInput(['time_splitting', 'use_state_mollifier'], fallback=False)):
            self.base = BaseCommanderExtended(config, self.scriptor, self.const, self.fl)
        else:
            self.base = BaseCommander(config, self.scriptor, self.const, self.fl)

        self.inputFile = config.filename
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

    def setupInitialSteps(self, zeroBaseline=True):
        if(self.const.saveStateLog):
            self.readStateDistance()

        commandString = ''
        J0, dummy = self.base.QoI()

        commands = []
        for k in range(1, self.const.Nsplit):
            commands += ['cp x0/%s ./a/ ' % self.fl.icFiles[k]]
        if (not zeroBaseline):
            for j in range(self.const.NcontrolRegion):
                commands += ['cp x0/%s ./a/ ' % self.const.globalControlSpaceFiles[j]]
        commandString += self.scriptor.nonMPILoopCommand(commands,'copy_control_params')

        numFiles = len(os.listdir(self.fl.LINMINDIR))
        if (numFiles>0):
            df_file = '%s/%d/%s' % (self.fl.LINMINDIR, numFiles-1, self.fl.lineMinLog)
            df_temp = pd.read_csv(df_file, sep='\t', header=0)
            lastStep = np.array(df_temp['step'][df_temp['directory index']=='b'])[0]

        steps, Js = np.zeros(2), np.zeros(2)
        steps[1] = self.const.initial_step if numFiles==0 else lastStep
        Js[0] = J0

        temp = ['x0/'+file for file in self.fl.globalControlSpaceFiles]
        if (zeroBaseline):
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

        command = 'python3 ' + self.fl.decisionMaker + ' ' + self.inputFile + ' 3'
        if (zeroBaseline):
            command += ' -zero_baseline'
        self.writeCommandFile(command+'\n', self.fl.decisionMakerCommandFile)

        data = {'step':steps,'QoI':Js,'directory index':['a','x']}
        df = pd.DataFrame(data)
        df.to_csv(self.fl.lineMinLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)
        print ('Initial steps written in command. Run ' + self.fl.globalCommandFile + '.')
        return 0

    def nextMnbrak(self, zeroBaseline=True):
        df = self.base.collectQoI(self.fl.lineMinLog)
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
            print ('Forward run crashed. Retrying with a smaller step.')
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

                self.writeCommandFile('exit 0\n', self.fl.globalCommandFile)
                command = 'python3 ' + self.fl.decisionMaker + ' ' + self.inputFile + ' 4 -linmin_initial'
                if(zeroBaseline):
                    command += ' -zero_baseline'
                self.writeCommandFile(command+'\n', self.fl.decisionMakerCommandFile)

                df.to_csv(self.fl.lineMinLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)
                print (np.array(df[df['directory index']!='0']))
                print ('MNBRAK: initial mininum bracket is prepared.')
                return 0
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
                print ('MNBRAK: expanding the bracket - Run intermediate forward simulation.')
        elif (cExists):
            c = np.array(df['step'][df['directory index']=='c'])[0]
            Jc = np.array(df['QoI'][df['directory index']=='c'])[0]

            if (Jx < Ja):
                self.base.switchDirectory('x','b',df)

                self.writeCommandFile('exit 0\n', self.fl.globalCommandFile)
                command = 'python3 ' + self.fl.decisionMaker + ' ' + self.inputFile + ' 4 -linmin_initial'
                if (zeroBaseline):
                    command += ' -zero_baseline'
                self.writeCommandFile(command+'\n', self.fl.decisionMakerCommandFile)

                df.to_csv(self.fl.lineMinLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)
                print (np.array(df[df['directory index']!='0']))
                print ('MNBRAK: initial mininum bracket is prepared.')
                return 0
            else:
                self.base.switchDirectory('x','c',df)
                df.loc[df['directory index']=='x','directory index'] = '0'

                new_x = x / self.const.golden_ratio
                print ('MNBRAK: narrowing the bracket - Run intermediate forward simulations.')
        else:
            if (Jx > Ja):
                self.base.switchDirectory('x','c',df)

                new_x = x / self.const.golden_ratio
                print ('MNBRAK: narrowing the bracket - Run intermediate forward simulations.')
            else:
                self.base.switchDirectory('x','b',df)

                if (len(df[(df['QoI'] > self.const.huge)]) > 0):
                    maxX = np.min(np.array(df['step'][df['QoI'] > self.const.huge]))
                    Cr = 1.0 - 1.0 / self.const.golden_ratio
                    new_x = x + Cr * (maxX-x)
                else:
                    new_x = x * self.const.golden_ratio if (x < self.const.safe_zone) else x * 1.05
                print ('MNBRAK: expanding the bracket - Run intermediate forward simulation.')

        data = {'step':[new_x], 'QoI':[np.nan], 'directory index':['x']}
        new_df = pd.DataFrame(data)
        df = df.append(new_df, ignore_index=True)

        commandString = ''
        temp = ['x0/'+file for file in self.fl.globalControlSpaceFiles]
        if (zeroBaseline):
            for j in range(self.const.NcontrolRegion):
                temp[j] = ''
        target = ['x/'+file for file in self.fl.globalControlSpaceFiles]
        commandString += self.scriptor.parallelPurgeCommand(target,'purge_target')
        commandString += self.base.zaxpyCommand(target, -new_x, self.fl.globalConjugateGradientFiles, temp)
        commandString += '\n'
        commandString += self.base.forwardRunCommand('x')
        commandString += self.scriptor.purgeDirectoryCommand('x')
        self.writeCommandFile(commandString, self.fl.globalCommandFile)

        command = 'python3 ' + self.fl.decisionMaker + ' ' + self.inputFile + ' 3'
        if (zeroBaseline):
            command += ' -zero_baseline'
        self.writeCommandFile(command+'\n', self.fl.decisionMakerCommandFile)

        df.to_csv(self.fl.lineMinLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)
        print (df[df['directory index']!='0'])

        return 0

    def nextLinmin(self, zeroBaseline=True, initial=True, stop=False):
        if (initial):
            df = pd.read_csv(self.fl.lineMinLog, sep='\t', header=0)
        else:
            df = self.base.collectQoI(self.fl.lineMinLog)

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
            df.to_csv(self.fl.lineMinLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)

        if (stop):
            print (df[df['directory index']!='0'])
            print ('LINMIN: line minimization is stopped.')

            self.writeCommandFile('exit 0\n', self.fl.globalCommandFile)
            command = 'python3 ' + self.fl.decisionMaker + ' ' + self.inputFile + ' 5'
            if(zeroBaseline):
                command += ' -zero_baseline'
            command += '\n exit 0\n'
            self.writeCommandFile(command, self.fl.decisionMakerCommandFile)
            return 0

        a = np.array(df['step'][df['directory index']=='a'])[0]
        b = np.array(df['step'][df['directory index']=='b'])[0]
        c = np.array(df['step'][df['directory index']=='c'])[0]
        Ja = np.array(df['QoI'][df['directory index']=='a'])[0]
        Jb = np.array(df['QoI'][df['directory index']=='b'])[0]
        Jc = np.array(df['QoI'][df['directory index']=='c'])[0]

        if ( (c-a) < b * self.const.linminTol + self.const.eps ):
            print ('steps: %.16E %.16E %.16E' % (a,b,c))
            print ('QoIs: %.16E %.16E %.16E' % (Ja,Jb,Jc))
            print ('LINMIN: line minimization is finished.')

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
            if( zeroBaseline and (not self.const.ignoreController) ):
                commandString += self.base.generalSetOptionCommand
                targetInputFiles = ['x0/%s/%s'%(dir,file) for dir, file in zip(self.fl.directories, self.fl.inputFiles)]
                commands = []
                for k in range(self.const.Nsplit):
                    commands += ['./setOption.sh %s "controller_switch" true' % targetInputFiles[k]]
                commandString += self.scriptor.nonMPILoopCommand(commands,'magudi_option_turn_on_controller')

            commandString += self.base.forwardRunCommand()

            self.writeCommandFile(commandString, self.fl.globalCommandFile)
            command = 'python3 ' + self.fl.decisionMaker + ' ' + self.inputFile + ' 5'
            if(zeroBaseline):
                command += ' -zero_baseline'                                        # afterLinmin doesn't use zeroBaseline. This has no effect.
            self.writeCommandFile(command+'\n', self.fl.decisionMakerCommandFile)
            return 0

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
        if (zeroBaseline):
            for j in range(self.const.NcontrolRegion):
                temp[j] = ''
        target = ['x/'+file for file in self.fl.globalControlSpaceFiles]
        commandString += self.scriptor.parallelPurgeCommand(target,'purge_target')
        commandString += self.base.zaxpyCommand(target, -xs, self.fl.globalConjugateGradientFiles, temp)
        commandString += '\n'
        commandString += self.base.forwardRunCommand('x')
        commandString += self.scriptor.purgeDirectoryCommand('x')
        self.writeCommandFile(commandString, self.fl.globalCommandFile)

        command = 'python3 ' + self.fl.decisionMaker + ' ' + self.inputFile + ' 4'
        if(zeroBaseline):
            command += ' -zero_baseline'
        self.writeCommandFile(command+'\n', self.fl.decisionMakerCommandFile)

        df.to_csv(self.fl.lineMinLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)
        print (df[df['directory index']!='0'])
        print ('LINMIN: next linmin evaluation is prepared-')
        return 1

    def beforeLinmin(self, initial, zeroBaseline):
        J0, subJ0 = self.base.QoI()
        gg, subgg = self.base.readInnerProduct(self.fl.ggFiles)

        data = [[J0]+list(subJ0)]
        J_new_df = pd.DataFrame(data,columns=self.fl.forwardLogColumns)
        data = [[gg]+list(subgg)]
        dJ_new_df = pd.DataFrame(data,columns=self.fl.gradientLogColumns)

        commandString = ''
        if (self.const.saveStateLog):
            commandString += self.base.checkStateDistance()

        if (initial):
            J_new_df.to_csv(self.fl.forwardLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)
            dJ_new_df.to_csv(self.fl.gradientLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)

            df = pd.DataFrame([np.nan],columns=['reduction'])
            df.to_csv(self.fl.CGLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)

            commands = []
            for k in range(self.const.NcontrolSpace):
                commands += ['cp %s %s' % (self.fl.globalGradFiles[k], self.fl.globalConjugateGradientFiles[k])]
            commandString += self.scriptor.nonMPILoopCommand(commands,'initial-conjugate-gradient')
            self.writeCommandFile(commandString, self.fl.globalCommandFile)

            command = 'python3 ' + self.fl.decisionMaker + ' ' + self.inputFile + ' 2'
            if (zeroBaseline):
                command += ' -zero_baseline'
            self.writeCommandFile(command+'\n', self.fl.decisionMakerCommandFile)

            print ('Initial line minimization is ready. Run mnbrak and linmin procedures.')
            return 0

        df = pd.read_csv(self.fl.forwardLog, sep='\t', header=0)
        J1 = df.at[df.index[-1],'total']
        df = df.append(J_new_df)
        df.to_csv(self.fl.forwardLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)

        df = pd.read_csv(self.fl.gradientLog, sep='\t', header=0)
        gg1 = df.at[df.index[-1],'total']
        gg0 = df.at[df.index[0],'total']
        df = df.append(dJ_new_df)
        df.to_csv(self.fl.gradientLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)

        df = pd.read_csv(self.fl.CGLog, sep='\t', header=0)
        df_addendum = pd.DataFrame([(1.0-J0/J1)], columns=['reduction'])
        df = df.append(df_addendum)
        df.to_csv(self.fl.CGLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)

        #reduction = gg/gg0
        #if (reduction<=tol):
        if (gg <= self.const.tol):
            print ('FRPRMN - after linmin: conjugate-gradient optimization is finished.')
            self.writeCommandFile('exit 1\n', self.fl.globalCommandFile)
            self.writeCommandFile('exit 1\n', self.fl.decisionMakerCommandFile)
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

        command = 'python3 ' + self.fl.decisionMaker + ' ' + self.inputFile + ' 2'
        if(zeroBaseline):
            command += ' -zero_baseline'
        self.writeCommandFile(command+'\n', self.decisionMakerCommandFile)

        print ('line minimization is ready. Run mnbrak and linmin procedures.')
        return 0

    def afterLinmin(self, zeroBaseline):
        #copy line minimization log
        numFiles = len(os.listdir(self.fl.LINMINDIR))
        subprocess.check_call('mkdir -p %s/%d' % (self.fl.LINMINDIR, numFiles), shell=True)
        subprocess.check_call('cp %s %s/%d/'%(self.fl.lineMinLog, self.fl.LINMINDIR, numFiles), shell=True)

        commandString = ''

        # copy result files
        commands = []
        for k in range(self.const.Nsplit):
            costFunctionalFile = 'b/%s/%s.cost_functional.txt' % (self.fl.directories[k], self.fl.prefixes[k])
            commands += ['cp %s linminLog/%d/' % (costFunctionalFile,numFiles)]
            commands += ['cp %s linminLog/%d/' % (self.fl.diffOutputFiles[k],numFiles)]

        commandString += self.scriptor.nonMPILoopCommand(commands,'saving_line_minimization_files')

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

        commands = []
        for k in range(self.const.Nsplit):
            costSensitivityFile = 'x0/%s/%s.cost_sensitivity.txt'%(self.fl.directories[k], self.fl.prefixes[k])
            commands += ['cp %s %s/%d/ '%(costSensitivityFile, self.fl.LINMINDIR, numFiles)]
        commandString += self.scriptor.nonMPILoopCommand(commands,'copy_gradient_history')

        # Polak-Ribiere
        commandString += self.base.dggCommand()
        commandString += self.scriptor.parallelPurgeCommand(self.fl.previousGradFiles,'purge_prev_grad')
        self.writeCommandFile(commandString, self.fl.globalCommandFile)
        command = 'python3 ' + self.fl.decisionMaker + ' ' + self.inputFile + ' 1 \n'
        if (self.const.useLagrangian and (numFiles > self.const.Nlinmin)):
            command += 'exit -1 \n'
        self.writeCommandFile(command, self.fl.decisionMakerCommandFile)

        print ('FRPRMN - after linmin: postprocessing is finished. Run new forward/adjoint simulations.')
        return 1
