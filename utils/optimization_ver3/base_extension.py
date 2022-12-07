from penalty_norm import penaltySwitcher

import numpy as np
import subprocess
import pandas as pd

from base import BaseCommander

__all__ = ['BaseCommanderExtended']

class BaseCommanderExtended(BaseCommander):

    def forwardRunCommand(self, baseDirectory='x0', zeroControlForcing=False):
        bdir = baseDirectory

        commandString = ''
        commandString += self.distributeCommand(bdir, zeroControlForcing)

        commandString += 'crash=0\n'

        commands, commandDirs = [], []
        for k in range(self.const.Nsplit):
            zFile = '%s/%s/%s' % (bdir, self.fl.directories[k], self.fl.icFiles[k])
            xFile = '%s/%s' % (bdir, self.fl.icFiles[k])
            if(k==0):
                commandString += 'cp %s %s \n' % (xFile, zFile)
            else:
                yFile = '%s/%s/%s' % (bdir, self.fl.directories[k-1], self.fl.matchingForwardFiles[k-1] )
                commands = ['./patchup_qfile %s %s %s %s' % (zFile, self.fl.icMollifierFile, xFile, yFile)]
                commandString += self.scriptor.singleJobCommand(commands,'qfile-zaxpy','patch_up_initial_condition%d' % k)

            commandDirs = ['%s/%s' % (bdir, self.fl.directories[k])]
            commands = ['./forward --input %s' % self.fl.inputFiles[k]]
            commandString += self.scriptor.singleJobCommand(commands,'forward','forward%d' % k, directories=commandDirs)

            commandDir = '%s/%s' % (bdir, self.fl.directories[k])
            crashCheck = 'if [ -f "%s/%s-crashed.q" ]; then\n' % (commandDir, self.fl.prefixes[k])
            crashCheck += '   echo "%d-th forward run crashed. Retry with a smaller step."\n' % k
            crashCheck += '   mv %s/%s-crashed.q %s/%s\n' % (commandDir, self.fl.prefixes[k],
                                                             commandDir, self.fl.matchingForwardFiles[k])
            crashCheck += '   let "crash+=1"\n'
            crashCheck += 'fi\n\n'
            commandString += crashCheck

        commandString += 'if [ $crash -eq 0 ]; then\n'                                              \
                          '   echo "No crash is found."\n'                                          \
                          'else\n'                                                                  \
                          '   echo "$crash crashes found."\n'                                       \
                          '   exit 0\n'                                                             \
                          'fi\n\n'

        if (self.const.terminalObjective):
            commandDirs = ['%s/%s' % (bdir, self.fl.directories[-1])]
            commands = ['./terminal_objective --mode forward --input %s --output %s'                \
                        % (self.fl.inputFiles[-1], self.fl.terminalOutputFile) ]
            commandString += self.scriptor.singleJobCommand(commands,'qfile_zaxpy',
                                                          'terminal_objective',directories=commandDirs)

        commands = []
        for k in range(self.const.Nsplit):
            matchingFile = '%s/%s/%s'%(bdir, self.fl.directories[k-1], self.fl.matchingForwardFiles[k-1])
            icFile = '%s/%s' % (bdir, self.fl.icFiles[k])
            commands += ['./qfile_zaxpy %s %.16E %s %s --input %s'                              \
                        % (self.fl.diffFiles[k], -1.0, matchingFile, icFile, self.fl.globalInputFile)]
        commandString += self.scriptor.parallelLoopCommand(commands,'qfile-zaxpy','forward_qfile_zaxpy')

        commands = []
        for k in range(self.const.Nsplit):
            commands += ['./spatial_inner_product %s %s --mollifier %s --output %s --input %s'                 \
                        % (self.fl.diffFiles[k], self.fl.diffFiles[k], self.fl.icMollifierFile,
                           self.fl.diffOutputFiles[k], self.fl.globalInputFile)]
        commandString += self.scriptor.parallelLoopCommand(commands,'qfile-zxdoty','forward_qfile_zxdoty')

        if (self.const.useLagrangian):
            commands = []
            for k in range(self.const.Nsplit):
                commands += ['./spatial_inner_product %s %s --mollifier %s --output %s --input %s'                 \
                            % (self.fl.diffFiles[k], self.fl.lagrangianFiles[k], self.fl.icMollifierFile,
                               self.fl.lagrangianOutputFiles[k], self.fl.globalInputFile)]
            commandString += self.scriptor.parallelLoopCommand(commands,'qfile-zxdoty','forward_qfile_lagrangian')

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

        commands = []
        for k in range(self.const.Nsplit):
            commands += ['./setOption.sh %s "adjoint_nonzero_initial_condition" "true"'                 \
                        % targetInputFiles[k]]
        commandString += self.scriptor.nonMPILoopCommand(commands,'magudi_option_nonzero_initial_condition_true')

        L2sqSum = 0.0
        kStart = 0 if self.const.periodicSolution else 1
        for k in range(kStart, self.const.Nsplit):
            L2sqSum += self.readScalar(self.fl.diffOutputFiles[k])

        commands = []
        for k in range(0,-self.const.Nsplit,-1):
            # Assuming all weights are equal except k=0.
            weight = -self.const.matchingConditionWeight[k] / self.penalty.gradientFactor(L2sqSum)
            matchingAdjointFile = '%s/%s/%s' % (bdir, self.fl.directories[k-1], self.fl.matchingAdjointFiles[k-1])
            commands = ['./qfile_zaxpy %s %.16E %s --zero --input %s'                                      \
            % (matchingAdjointFile, weight, self.fl.diffFiles[k], self.fl.globalInputFile)]
            commandString += self.scriptor.singleJobCommand(commands,'qfile-zaxpy','adjoint_run_qfile_zaxpy%d' % k)

            if (self.const.useLagrangian):
                matchingAdjointFile = '%s/%s/%s' % (bdir, self.fl.directories[k-1], self.fl.matchingAdjointFiles[k-1])
                commands = ['./qfile_zaxpy %s %.16E %s %s --input %s'                                      \
                % (matchingAdjointFile, 1.0, self.fl.lagrangianFiles[k], self.fl.matchingAdjointFile, self.fl.globalInputFile)]
                commandString += self.scriptor.singleJobCommand(commands,'qfile-zaxpy','adjoint_run_qfile_lagrangian%d' % k)

            matchingAdjointFile = '%s/%s/%s' % (bdir, self.fl.directories[k-1], self.fl.matchingAdjointFiles[k-1])
            if(k==0):
                patchingAdjointFile = '--zero'
            else:
                patchingAdjointFile = '%s/%s/%s' % (bdir, self.fl.directories[k], self.fl.icAdjointFiles[k])
            commands = ['./patchup_qfile %s %s %s %s --input %s'                                      \
            % (matchingAdjointFile, self.fl.icMollifierFile, matchingAdjointFile, patchingAdjointFile, self.fl.globalInputFile)]
            commandString += self.scriptor.singleJobCommand(commands,'qfile-zaxpy','adjoint_run_patch_up_qfile%d' % k)

            if ((k==0) and self.const.terminalObjective):
                commandDirs = ['%s/%s' % (bdir, self.fl.directories[-1])]
                commands = ['./terminal_objective --mode adjoint --input %s' % self.fl.inputFiles[-1]]
                commandString += self.scriptor.singleJobCommand(commands,'qfile_zaxpy',
                                                                'terminal_sensitivity',directories=commandDirs)

            commandDirs = ['%s/%s' % (bdir, self.fl.directories[k-1])]
            commands = ['./adjoint --input %s' % self.fl.inputFiles[k-1]]
            commandString += self.scriptor.singleJobCommand(commands,'adjoint',
                                                            'adjoint%d' % k,directories=commandDirs)

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
                commands += ['./qfile_zaxpy %s %.16E %s %s --input %s' % (self.fl.globalGradFiles[idx],
                            -1.0, self.fl.lagrangianFiles[k], self.fl.globalGradFiles[idx], self.fl.globalInputFile)]
            commandString += self.scriptor.parallelLoopCommand(commands,'qfile-zaxpy','adjoint_run_ic_lagrangian')

        commands = []
        for k in range(self.const.Nsplit):
            idx = self.const.NcontrolRegion + k
            commands += ['./patchup_qfile %s %s %s --zero --input %s' % (self.fl.globalGradFiles[idx],
                        self.fl.icMollifierFile, self.fl.globalGradFiles[idx], self.fl.globalInputFile)]
        commandString += self.scriptor.parallelLoopCommand(commands,'qfile-zaxpy','patch_up_adjoint_run_ic')

        return commandString
