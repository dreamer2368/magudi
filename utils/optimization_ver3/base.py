from penalty_norm import penaltySwitcher

import numpy as np
import subprocess
import pandas as pd

class BaseCommander:
    scriptor = None
    penalty = None
    const = None
    fl = None

    magudiSetOptionCommand = 'function setOption() {\n'                                         \
                             '    if grep -q "$1" magudi.inp\n'                                 \
                             '    then\n'                                                       \
                             '	sed -i "s+^.*$1.*$+$1 = $2+g" magudi.inp\n'                     \
                             '    else\n'                                                       \
                             '	echo "$1 = $2" >> magudi.inp\n'                                 \
                             '    fi\n'                                                         \
                             '}\n'

    generalSetOptionCommand = 'cat <<EOF > setOption.sh\n'                                   \
                              '#!/bin/bash\n'                                                \
                              'if grep -q "\$2" \$1 \n'                                      \
                              'then\n'                                                       \
                              'sed -i.bu "s+^.*\$2.*\$+\$2 = \$3+g" \$1 \n'                  \
                              'else\n'                                                       \
                              'echo "\$2 = \$3" >> \$1 \n'                                   \
                              'fi\n'                                                         \
                              'EOF\n'                                                        \
                              'chmod u+x setOption.sh\n'

    def __init__(self, config, scriptor, const, filenamelist):
        self.scriptor = scriptor
        self.const = const
        self.fl = filenamelist

        self.penalty = penaltySwitcher.get(self.const.penaltyType)
        return

    def readScalar(self, scalarFilename):
        try:
            fID = open(scalarFilename,'r')
        except FileNotFoundError:
            print (scalarFilename+' is not found.')
            return np.nan
        else:
            scalar = float(fID.read())
            fID.close()
            return scalar

    def QoI(self, baseDirectory = 'x0'):
        bdir = baseDirectory

        subJ = np.zeros(3 * self.const.Nsplit + 1)
        J = 0.0
        for k in range(self.const.Nsplit):
            subJ[k] = readScalar(bdir + '/' + self.fl.directories[k] + '/' + self.fl.outputFiles[k])
            if (not self.const.ignoreIntegralObjective):
                J += subJ[k]

            if (self.const.useLagrangian):
                subJ[2 * self.const.Nsplit + k] = readScalar(self.fl.lagrangianOutputFiles[k])
                # followed the convention of Portryagin's minimum principle
                J -= subJ[2 * self.const.Nsplit + k]

        if (self.const.terminalObjective):
            subJ[3 * self.const.Nsplit] = readScalar(self.fl.terminalOutputFile)
            J += subJ[3 * self.const.Nsplit]

        L2sqSum = 0.0
        kStart = 0 if self.const.periodicSolution else 1
        for k in range(kStart, self.const.Nsplit):
            L2sq = readScalar(self.fl.diffOutputFiles[k])
            subJ[self.const.Nsplit + k] = 0.5 * L2sq
            L2sqSum += L2sq
        # Assuming all weights are equal except k=0.
        # J += self.const.matchingConditionWeight[-1] * self.penalty.norm(L2sqSum)
        J += np.dot(self.const.matchingConditionWeight, self.penalty.norm(L2sqSum))

        return J, subJ

    def collectQoI(self, logFilename):
        import pandas as pd

        new_J, dummy = self.QoI('x')

        df = pd.read_csv(logFilename, sep='\t', header=0)
        df.at[df['directory index']=='x','QoI'] = new_J
        df = df.sort_values(by='step',ascending=True)
        df = df.reset_index(drop=True)
        return df

    def innerProductCommand(self, xFiles, yFiles, outputFiles_=None):
        if (outputFiles_ is None):
            outputFiles_ = self.fl.innerProductFiles

        Nx, Ny = len(xFiles), len(yFiles)
        if ((Nx != self.const.NcontrolSpace) or (Ny != self.const.NcontrolSpace)):
            raise LookupError('Number of files must be equal to '+str(self.const.NcontrolSpace)+'!')

        commandString = ''

        commands = []
        for j in range(self.const.NcontrolRegion):
            commands += ['./zxdoty %s %s %s %s'                                                           \
                        % (outputFiles_[j], xFiles[j], yFiles[j], self.fl.globalNormFiles[j])]
        commandString += self.scriptor.parallelLoopCommand(commands,'zxdoty','zxdoty_control_forcing')

        commands = []
        for k in range(self.const.NcontrolRegion, self.const.NcontrolSpace):
            index = k - self.const.NcontrolRegion
            commands += ['./spatial_inner_product %s %s --input %s --output %s'                    \
                        % (xFiles[k], yFiles[k], self.fl.globalInputFile, outputFiles_[k])]
        commandString += self.scriptor.parallelLoopCommand(commands,'qfile-zxdoty',
                                                          'zxdoty_initial_condition')
        return commandString

    def readInnerProduct(self, outputFiles_=None):
        if (outputFiles_ is None):
            outputFiles_ = self.fl.innerProductFiles

        scalar = 0.0
        subScalar = np.zeros(self.const.NcontrolSpace)
        for k in range(self.const.NcontrolSpace):
            subScalar[k] = self.readScalar(outputFiles_[k])
        subScalar[-self.const.Nsplit:] *= self.const.initialConditionControllability
        scalar = np.sum(subScalar)

        return scalar, subScalar

    def zaxpyCommand(self, zFiles, a, xFiles, yFiles=None):
        Nx, Nz = len(xFiles), len(zFiles)
        if ((Nx != self.const.NcontrolSpace) or (Nz != self.const.NcontrolSpace)):
            raise LookupError('Number of files must be equal to '+str(self.const.NcontrolSpace)+'!')
        if (yFiles is not None):
            Ny = len(yFiles)
            if (Ny != self.const.NcontrolSpace):
                raise LookupError('Number of files must be equal to '+str(self.const.NcontrolSpace)+'!')

        commandString = ''

        commands = []
        for j in range(self.const.NcontrolRegion):
            command = './zaxpy ' + zFiles[j] + ' ' + "{:.16E}".format(a) + ' ' + xFiles[j]
            if (yFiles is not None):
                command += ' '+yFiles[j]
            commands += [command]
        commandString += self.scriptor.parallelLoopCommand(commands,'zaxpy',
                                                          'zaxpy_control_forcing')

        commands = []
        for k in range(self.const.NcontrolRegion, self.const.NcontrolSpace):
            index = k - self.const.NcontrolRegion
            w2 = self.const.initialConditionControllability[index]
            command = './qfile_zaxpy ' + zFiles[k] + ' ' + "{:.16E}".format(a * w2) + ' ' + xFiles[k]
            if (yFiles is not None):
                command += ' '+yFiles[k]
            commands += [command]
        commandString += self.scriptor.parallelLoopCommand(commands,'qfile-zaxpy',
                                                          'zaxpy_initial_condition')
        return commandString

    def distributeCommand(self, baseDirectory, zeroControlForcing):
        bdir = baseDirectory

        commandString = ''

        if (not zeroControlForcing):
            commands = []
            if (self.const.Nsplit==1):
                for j in range(self.const.NcontrolRegion):
                    sliceControlForcingFile = '%s/%s/%s%s' % (bdir, self.fl.directories[0],
                                                              self.fl.prefixes[0], self.fl.controlForcingFiles[j])
                    globalControlForcingFile = '%s/%s' % (bdir, self.fl.globalControlSpaceFiles[j])
                    commands += ['./zaxpy %s 1.0 %s' % (sliceControlForcingFile, globalControlForcingFile)]
                commandString += self.scriptor.parallelLoopCommand(commands,'zaxpy','slice_control_forcing')
            else:
                for k in range(self.const.Nsplit):
                    kOffset = self.const.Nts * k
                    for j in range(self.const.NcontrolRegion):
                        sliceControlForcingFile = '%s/%s/%s%s' % (bdir, self.fl.directories[k],
                                                                  self.fl.prefixes[k], self.fl.controlForcingFiles[j])
                        globalControlForcingFile = '%s/%s' % (bdir, self.fl.globalControlSpaceFiles[j])
                        commands += ['./slice_control_forcing %s %s %d %d %d'                                                           \
                                    % (sliceControlForcingFile, globalControlForcingFile, self.const.totalTimestep, kOffset, self.const.Nts)]
                commandString += self.scriptor.parallelLoopCommand(commands,'slice','slice_control_forcing')

        return commandString

    def gatherControlForcingGradientCommand(self):
        commandString = ''

        from os import path
        for j in range(self.const.NcontrolRegion):
            if (path.exists(self.fl.globalGradFiles[j])):
                print ("Previous global gradient file %s still exists. Purging it for safety." % (self.fl.globalGradFiles[j]))
                commandString += 'rm %s \n' % (self.fl.globalGradFiles[j])

        commands = []
        if (self.const.Nsplit==1):
            for j in range(self.const.NcontrolRegion):
                sliceGradFile = 'x0/%s/%s%s' % (self.fl.directories[0], self.fl.prefixes[0], self.fl.gradientFiles[j])
                commands += ['mv %s %s' % (sliceGradFile, self.fl.globalGradFiles[j])]
            commandString += self.scriptor.nonMPILoopCommand(commands,'paste_control_forcing')
        else:
            for k in range(self.const.Nsplit):
                kOffset = self.const.Nts * k
                for j in range(self.const.NcontrolRegion):
                    sliceGradFile = 'x0/%s/%s%s' % (self.fl.directories[k], self.fl.prefixes[k], self.fl.gradientFiles[j])
                    commands += ['./paste_control_forcing %s %s %d %d %d'                                                           \
                                % (self.fl.globalGradFiles[j], sliceGradFile, self.const.totalTimestep, kOffset, self.const.Nts)]
            commandString += self.scriptor.parallelLoopCommand(commands,'paste','paste_control_forcing')

        return commandString

    def switchDirectory(self, firstDirectory, secondDirectory, df=None):
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

    def forwardRunCommand(self, baseDirectory='x0', zeroControlForcing=False):
        bdir = baseDirectory

        commandString = ''
        commandString += self.distributeCommand(bdir, zeroControlForcing)


        commands, commandDirs = [], []
        for k in range(self.const.Nsplit):
            commandDirs += ['%s/%s' % (bdir, self.fl.directories[k])]
            commands += ['./forward --input %s' % self.fl.inputFiles[k]]
        commandString += self.scriptor.parallelLoopCommand(commands,'forward',
                                                          'forward',directories=commandDirs)

        crashCheck  = 'crash=0\n'
        for k in range(self.const.Nsplit):
            commandDir = '%s/%s' % (bdir, self.fl.directories[k])
            crashCheck += 'if [ -f "%s/%s-crashed.q" ]; then\n' % (commandDir, self.fl.prefixes[k])
            crashCheck += '   echo "%d-th forward run crashed. Retry with a smaller step."\n' % k
            crashCheck += '   mv %s/%s-crashed.q %s/%s\n' % (commandDir, self.fl.prefixes[k], commandDir, self.fl.matchingForwardFiles[k])
            crashCheck += '   let "crash+=1"\n'
            crashCheck += 'fi\n\n'
        crashCheck += 'if [ $crash -eq 0 ]; then\n'                                             \
                      '   echo "No crash is found."\n'                                          \
                      'else\n'                                                                  \
                      '   echo "$crash crashes found."\n'                                       \
                      '   exit 0\n'                                                             \
                      'fi\n\n'
        commandString += crashCheck

        if (self.const.terminalObjective):
            commandDirs = ['%s/%s' % (bdir, self.fl.directories[-1])]
            commands = ['./terminal_objective --mode forward --input %s --output %s'                \
                        % (self.fl.inputFiles[-1], self.fl.terminalOutputFile) ]
            commandString += self.scriptor.singleJobCommand(commands,'qfile-zaxpy',
                                                          'terminal_objective',directories=commandDirs)

        commands = []
        for k in range(self.const.Nsplit):
            matchingFile = '%s/%s/%s' % (bdir, self.fl.directories[k-1], self.fl.matchingForwardFiles[k-1])
            icFile = '%s/%s' % (bdir, self.fl.icFiles[k])
            commands += ['./qfile_zaxpy %s %.16E %s %s --input %s'                              \
                        % (self.fl.diffFiles[k], -1.0, matchingFile, icFile, self.fl.globalInputFile)]
        commandString += self.scriptor.parallelLoopCommand(commands,'qfile-zaxpy',
                                                          'forward_qfile_zaxpy')

        commands = []
        for k in range(self.const.Nsplit):
            commands += ['./spatial_inner_product %s %s --output %s --input %s'                 \
                        % (self.fl.diffFiles[k], self.fl.diffFiles[k], self.fl.diffOutputFiles[k], self.fl.globalInputFile)]
        commandString += self.scriptor.parallelLoopCommand(commands,'qfile-zxdoty',
                                                          'forward_qfile_zxdoty')

        if (self.const.useLagrangian):
            commands = []
            for k in range(self.const.Nsplit):
                commands += ['./spatial_inner_product %s %s --output %s --input %s'                 \
                            % (self.fl.diffFiles[k], self.fl.lagrangianFiles[k], self.fl.lagrangianOutputFiles[k], self.fl.globalInputFile)]
            commandString += self.scriptor.parallelLoopCommand(commands,'qfile-zxdoty',
                                                          'forward_qfile_lagrangian')

        return commandString

    def adjointRunCommand(self, baseDirectory='x0'):
        bdir = baseDirectory

        commandString = ''
        commandString += self.generalSetOptionCommand

        targetInputFiles = ['%s/%s/%s'%(bdir,dir,file) for dir, file in zip(self.fl.directories, self.fl.inputFiles)]

        if (self.const.ignoreIntegralObjective):
            commands = []
            for k in range(self.const.Nsplit):
                commands += ['./setOption.sh %s "adjoint_forcing_switch" "false"' % targetInputFiles[k]]
            commandString += self.scriptor.nonMPILoopCommand(commands,'magudi_option_disable_adjoint_forcing')

        L2sqSum = 0.0
        kStart = 0 if self.const.periodicSolution else 1
        for k in range(kStart, self.const.Nsplit):
            L2sqSum += readScalar(self.fl.diffOutputFiles[k])

        commands = []
        for k in range(self.const.Nsplit):
            # Assuming all weights are equal except k=0.
            weight = -self.const.matchingConditionWeight[k] / self.penalty.gradientFactor(L2sqSum)
            matchingAdjointFile = '%s/%s/%s' % (bdir, self.fl.directories[k-1], self.fl.matchingAdjointFiles[k-1])
            commands += ['./qfile_zaxpy %s %.16E %s --zero --input %s'                                      \
            % (matchingAdjointFile, weight, self.fl.diffFiles[k], self.fl.globalInputFile)]
        commandString += self.scriptor.parallelLoopCommand(commands,'qfile-zaxpy',
                                                            'adjoint_run_qfile_zaxpy')

        if (self.const.useLagrangian):
            commands = []
            for k in range(self.const.Nsplit):
                matchingAdjointFile = '%s/%s/%s' % (bdir, self.fl.directories[k-1], self.fl.matchingAdjointFiles[k-1])
                commands += ['./qfile_zaxpy %s %.16E %s %s --input %s'                                      \
                % (self.fl.matchingAdjointFile, 1.0, self.fl.lagrangianFiles[k], matchingAdjointFile, self.fl.globalInputFile)]
            commandString += self.scriptor.parallelLoopCommand(commands,'qfile-zaxpy',
                                                                'adjoint_run_qfile_lagrangian')

        commands = []
        for k in range(self.const.Nsplit):
            commands += ['./setOption.sh %s "adjoint_nonzero_initial_condition" "true"'                 \
                        % targetInputFiles[k]]
        commandString += self.scriptor.nonMPILoopCommand(commands,'magudi_option_nonzero_initial_condition_true')

        if (self.const.terminalObjective):
            commandDirs = ['%s/%s' % (bdir, self.fl.directories[-1])]
            commands = ['./terminal_objective --mode adjoint --input %s' % self.fl.inputFiles[-1]]
            commandString += self.scriptor.singleJobCommand(commands,'qfile-zaxpy',
                                                            'terminal_sensitivity',directories=commandDirs)

        commands, commandDirs = [], []
        for k in range(self.const.Nsplit):
            commandDirs += ['%s/%s' % (bdir, self.fl.directories[k])]
            commands += ['./adjoint --input %s' % self.fl.inputFiles[k]]
        commandString += self.scriptor.parallelLoopCommand(commands,'adjoint',
                                                            'adjoint',directories=commandDirs)

        # TODO: Either check the last adjoint run, or complete the periodic optimization!!
        if (self.const.ignoreIntegralObjective):
            commands = []
            for k in range(self.const.Nsplit):
                commands += ['./setOption.sh %s "adjoint_forcing_switch" "true"' % targetInputFiles[k]]
            commandString += self.scriptor.nonMPILoopCommand(commands,'magudi_option_enable_adjoint_forcing')

        commands = []
        for k in range(self.const.Nsplit):
            commands += ['./setOption.sh %s "adjoint_nonzero_initial_condition" "false"'                            \
                        % targetInputFiles[k]]
        commandString += self.scriptor.nonMPILoopCommand(commands,'magudi_option_nonzero_initial_condition_false')

        commands = []
        for k in range(self.const.Nsplit):
            idx = self.const.NcontrolRegion + k
            # Assuming all weights are equal except k=0.
            weight = self.const.matchingConditionWeight[k] / self.penalty.gradientFactor(L2sqSum)
            icAdjointFile = '%s/%s/%s' % (bdir, self.fl.directories[k], self.fl.icAdjointFiles[k])
            commands += ['./qfile_zaxpy %s %.16E %s %s --input %s'                                      \
                         % (self.fl.globalGradFiles[idx], weight, self.fl.diffFiles[k], icAdjointFile, self.fl.globalInputFile)]
        commandString += self.scriptor.parallelLoopCommand(commands,'qfile-zaxpy','adjoint_run_ic')

        if (self.const.useLagrangian):
            commands = []
            for k in range(self.const.Nsplit):
                idx = self.const.NcontrolRegion + k
                commands += ['./qfile_zaxpy %s %.16E %s %s --input %s'                                     \
                             % (self.fl.globalGradFiles[idx], -1.0, self.fl.lagrangianFiles[k], self.fl.globalGradFiles[idx], self.fl.globalInputFile)]
            commandString += self.scriptor.parallelLoopCommand(commands,'qfile-zaxpy','adjoint_run_ic_lagrangian')

        return commandString

    def dggCommand(self):
        commandString = ''

        commands = []
        for j in range(self.const.NcontrolRegion):
            commands += ['./zwxmwy %s %s %s %s %s'                                                    \
            %(self.fl.dggFiles[j], self.fl.globalGradFiles[j], self.fl.globalGradFiles[j], self.fl.previousGradFiles[j], self.fl.globalNormFiles[j])]
        commandString += self.scriptor.parallelLoopCommand(commands,'zxdoty','dgg_zwxmwy')

        commands = []
        for k in range(self.const.Nsplit):
            idx = self.const.NcontrolRegion + k
            diffFile = '%s/%s.diff.adjoint.q' % (self.fl.DIFFDIR, self.fl.prefixes[k])
            prevICGradFile = '%s' % (self.fl.previousGradFiles[idx])
            icGradFile = '%s' % (self.fl.globalGradFiles[idx])
            commands += ['./qfile_zaxpy %s %.16E %s %s --input %s'                                    \
                        % (diffFile, -1.0, prevICGradFile, icGradFile, self.fl.globalInputFile)]
        commandString += self.scriptor.parallelLoopCommand(commands,'qfile-zaxpy','dgg_qfile_zaxpy')

        commands = []
        for k in range(self.const.Nsplit):
            idx = self.const.NcontrolRegion + k
            diffFile = '%s/%s.diff.adjoint.q' % (self.fl.DIFFDIR, self.fl.prefixes[k])
            commands += ['./spatial_inner_product %s %s --output %s --input %s'                                \
                        % (diffFile, self.fl.globalGradFiles[idx], self.fl.dggFiles[idx], self.fl.globalInputFile)]
        commandString += self.scriptor.parallelLoopCommand(commands,'qfile-zxdoty','dgg_qfile_zxdoty')

        return commandString

    def updateLagrangian(self,weight,initial=False):
        commandString = ''

        commands = []
        for k in range(self.const.Nsplit):
            a = weight[k]
            command = './qfile_zaxpy ' + self.fl.lagrangianFiles[k] + ' ' + "{:.16E}".format(-a) + ' ' + self.fl.diffFiles[k]
            if (initial):
                command += ' --zero'
            else:
                command += ' '+self.fl.lagrangianFiles[k]
            commands += [command]
        commandString += self.scriptor.parallelLoopCommand(commands,'qfile-zaxpy','update_lagrangian')

        return commandString

    def checkStateDistance(self):
        commandString = ''

        commands = []
        for k in range(self.const.Nsplit):
            icFile = 'x0/%s' % self.fl.icFiles[k]
            baselineFile = self.fl.baselineFiles[k]
            baseDiffFile = '%s/%s.base_diff.q' % (self.fl.DIFFDIR, self.fl.prefixes[k])
            commands += ['./qfile_zaxpy %s -1.0 %s %s' % (baseDiffFile, baselineFile, icFile)]
        commandString += self.scriptor.parallelLoopCommand(commands,'qfile-zaxpy','baseline_diff')

        commands = []
        for k in range(self.const.Nsplit):
            diffFile = self.fl.diffFiles[k]
            baseDiffFile = '%s/%s.base_diff.q' % (self.fl.DIFFDIR, self.fl.prefixes[k])
            for j in range(NdiffRegion):
                diffOutput = '%s/%s.diff.%d.txt' % (self.fl.TXTDIR, self.fl.prefixes[k], j)
                baseDiffOutput = '%s/%s.base_diff.%d.txt' % (self.fl.TXTDIR, self.fl.prefixes[k], j)
                mollifierFile = self.fl.diffMollifierFiles[j]
                commands += ['./spatial_inner_product %s %s --mollifier %s --output %s --input %s'                 \
                            % (diffFile, diffFile, mollifierFile, diffOutput, self.fl.globalInputFile)]
                commands += ['./spatial_inner_product %s %s --mollifier %s --output %s --input %s'                 \
                            % (baseDiffFile, baseDiffFile, mollifierFile, baseDiffOutput, self.fl.globalInputFile)]
        commandString += self.scriptor.parallelLoopCommand(commands,'qfile-zxdoty','mollified_diffs')

        return commandString

    def readStateDistance(self):
        diffs = np.zeros(2 * self.const.Nsplit * self.const.NdiffRegion)
        for k in range(self.const.Nsplit):
            for j in range(self.const.NdiffRegion):
                diffOutput = '%s/%s.diff.%d.txt' % (self.fl.TXTDIR, self.fl.prefixes[k],j)
                baseDiffOutput = '%s/%s.base_diff.%d.txt' % (self.fl.TXTDIR, self.fl.prefixes[k], j)
                idx = 2 * (j + NdiffRegion * k)
                diffs[idx + 1] = self.readScalar(diffOutput)
                diffs[idx] = self.readScalar(baseDiffOutput)

        new_df = pd.DataFrame([list(diffs)],columns=self.fl.stateLogColumns)
        from os import path
        if (path.exists(self.fl.stateLog)):
            df = pd.read_csv(self.fl.stateLog, sep='\t', header=0)
            df = df.append(new_df)
        else:
            df = new_df
        df.to_csv(self.fl.stateLog, float_format='%.16E', encoding='utf-8', sep='\t', mode='w', index=False)

        return
