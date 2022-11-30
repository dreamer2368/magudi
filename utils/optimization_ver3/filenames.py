# from constants import *

class FilenameList:
    const = None

    decisionMaker = ''
    globalInputFile = ''

    globalPrefix = ''
    globalCommandFile = ''
    decisionMakerCommandFile = ''
    lineMinLog = ''
    CGLog = ''

    ROOTDIR = ''
    OUTDIR = ''
    DIFFDIR = ''
    GRADDIR = ''
    LGRNGDIR = ''
    TXTDIR = ''
    CGDIR = ''
    PREVDIR = ''
    DIFFLOGDIR = ''
    BASELINEDIR = ''
    LINMINDIR = ''
    dirList = ['']

    forwardLog = ''
    forwardLogColumns = ['']

    gradientLog = ''
    gradientLogColumns = ['']

    stateLog = ''
    stateLogColumns = ['']

    controlRegions = ['']

    gradientFiles = ['']
    normFiles = ['']
    conjugateGradientFiles = ['']
    controlForcingFiles = ['']

    terminalOutputFile = ''
    icMollifierFile = ''
    diffMollifierFiles = ['']

    prefixes = ['']
    inputFiles = ['']
    outputFiles = ['']
    commandFiles = ['']
    icFiles = ['']
    diffFiles = ['']
    diffOutputFiles = ['']
    matchingForwardFiles = ['']
    matchingAdjointFiles = ['']
    icAdjointFiles = ['']
    icGradientFiles = ['']
    directories = ['']
    lagrangianFiles = ['']
    lagrangianOutputFiles = ['']
    baselineFiles = ['']

    innerProductFiles = ['']
    ggFiles = ['']
    dggFiles = ['']

    globalGradFiles = ['']
    globalNormFiles = ['']
    globalConjugateGradientFiles = ['']
    globalControlSpaceFiles = ['']
    globalSensitivityFiles = ['']

    previousGradFiles = ['']
    previousCGFiles = ['']

    def __init__(self, config, const):
        self.const = const
        self.decisionMaker = 'optimization.py'

        self.globalInputFile = 'magudi.inp'

        ##################           optimization files              #####################

        self.globalPrefix = config.getInput(['global_prefix'], datatype=str)
        self.globalCommandFile = '%s.command.sh' % self.globalPrefix
        self.decisionMakerCommandFile = '%s.command.python.ready.sh' % self.globalPrefix
        self.lineMinLog = self.globalPrefix + '.line_minimization.txt'
        self.CGLog = self.globalPrefix + '.conjugate_gradient.txt'

        from os import path, getcwd
        self.ROOTDIR = getcwd()
        self.OUTDIR = '%s/out' % self.ROOTDIR
        self.DIFFDIR = '%s/diff' % self.ROOTDIR
        self.GRADDIR = '%s/grad' % self.ROOTDIR
        self.LGRNGDIR = '%s/lagrangian' % self.ROOTDIR
        self.TXTDIR = '%s/txt' % self.ROOTDIR
        self.CGDIR = '%s/cg' % self.ROOTDIR
        self.PREVDIR = '%s/previous' % self.ROOTDIR
        self.DIFFLOGDIR = '%s/diffLog' % self.ROOTDIR
        self.BASELINEDIR = '%s/baseline' % self.ROOTDIR
        self.LINMINDIR = '%s/linminLog' % self.ROOTDIR
        self.dirList = [self.OUTDIR, self.DIFFDIR, self.GRADDIR, self.LGRNGDIR,
                        self.TXTDIR, self.CGDIR, self.PREVDIR, self.DIFFLOGDIR,
                        self.LINMINDIR]

        self.forwardLog = self.globalPrefix + '.cost_reduction.txt'
        self.forwardLogColumns = ['total']
        for k in range(self.const.Nsplit):
            self.forwardLogColumns += ['cost_functional%d' % k]
        for k in range(self.const.Nsplit):
            self.forwardLogColumns += ['matching_penalty%d' % k]
        for k in range(self.const.Nsplit):
            self.forwardLogColumns += ['lagrangian%d' % k]
        self.forwardLogColumns += ['terminal']

        self.gradientLog = self.globalPrefix + '.gg_log.txt'
        self.gradientLogColumns = ['total']
        for j in range(self.const.NcontrolRegion):
            self.gradientLogColumns += ['control_forcing%d' % j]
        for k in range(self.const.Nsplit):
            self.gradientLogColumns += ['initial_condition%d' % k]

        self.stateLog = self.globalPrefix + '.state_log.txt'
        self.stateLogColumns = []
        for k in range(self.const.Nsplit):
            for j in range(self.const.NdiffRegion):
                self.stateLogColumns += ['base_diff%d_%d' % (k, j)]
                self.stateLogColumns += ['diff%d_%d' % (k, j)]

        ##################           global level files              #####################

        self.controlRegions = config.getInput(['controller', 'actuator_list'], fallback=[])
        assert(len(self.controlRegions) == self.const.NcontrolRegion)

        self.gradientFiles, self.normFiles, self.conjugateGradientFiles, self.controlForcingFiles = [], [], [], []
        for k in range(self.const.NcontrolRegion):
            self.gradientFiles += ['.gradient_%s.dat' % self.controlRegions[k]]
            self.normFiles += ['.norm_%s.dat' % self.controlRegions[k]]
            self.conjugateGradientFiles += ['.conjugate_gradient_%s.dat' % self.controlRegions[k]]
            self.controlForcingFiles += ['.control_forcing_%s.dat' % self.controlRegions[k]]

        self.terminalOutputFile = self.TXTDIR + '/' + self.globalPrefix + '.terminal_objective.txt'

        self.icMollifierFile = self.ROOTDIR + '/' + self.globalPrefix + '.ic_mollifier.f'
        self.diffMollifierFiles = [self.ROOTDIR + '/' + self.globalPrefix + '.ic_mollifier.%d.f' % k for k in range(self.const.NdiffRegion)]

        ##################           slice level files              #####################

        self.prefixes, self.inputFiles, self.outputFiles, self.commandFiles = [], [], [], []
        self.icFiles, self.diffFiles, self.diffOutputFiles                  = [], [], []
        self.matchingForwardFiles, self.matchingAdjointFiles                = [], []
        self.icAdjointFiles, self.icGradientFiles, self.directories         = [], [], []
        self.lagrangianFiles, self.lagrangianOutputFiles                    = [], []
        self.baselineFiles                                                  = []
        for k in range(self.const.Nsplit):
            self.prefixes += ['%s-%01d' % (self.globalPrefix, k)]
            self.inputFiles += ['magudi-%01d.inp' % (k)]
            self.outputFiles += ['%s-%01d.forward_run.txt' % (self.globalPrefix, k)]
            self.commandFiles += ['%s-%01d.command.sh' % (self.globalPrefix, k)]
            self.icFiles += ['%s-%01d.ic.q' % (self.globalPrefix, k)]

            kOffset = self.const.startTimestep + self.const.Nts * (k + 1)
            self.diffFiles += ['%s-%01d.diff.q' % (self.globalPrefix, k)]                              #use only 0, ... , Nsplit-2
            self.diffOutputFiles += ['%s-%01d.diff.txt' % (self.globalPrefix, k)]                      #use only 0, ... , Nsplit-2
            self.matchingForwardFiles += ['%s-%01d-%08d.q' % (self.globalPrefix, k, kOffset)]          #use only 0, ... , Nsplit-2
            self.matchingAdjointFiles += ['%s-%01d-%08d.adjoint.q' % (self.globalPrefix, k, kOffset)]  #use only 0, ... , Nsplit-2
            self.lagrangianFiles += ['%s-%01d.lagrangian_multiplier.q' % (self.globalPrefix, k)]       #use only 0, ... , Nsplit-2
            self.lagrangianOutputFiles += ['%s-%01d.lagrangian.txt' % (self.globalPrefix, k)]          #use only 0, ... , Nsplit-2

            kOffset = self.const.startTimestep + self.const.Nts * k
            self.icAdjointFiles += ['%s-%01d-%08d.adjoint.q' % (self.globalPrefix, k, kOffset)]        #use only 1, ... , Nsplit-1
            self.icGradientFiles += ['%s-%01d.ic.adjoint.q' % (self.globalPrefix, k)]                  #use only 1, ... , Nsplit-1
            self.baselineFiles += ['%s-%08d.q' % (self.globalPrefix, kOffset)]
            self.directories += ['%01d' % k]

        ##################           global+slice level files              #####################

        self.innerProductFiles, self.ggFiles, self.dggFiles = [], [], []
        for j in range(self.const.NcontrolRegion):
            self.innerProductFiles += ['%s.inner_product_%s.txt' % (self.globalPrefix, self.controlRegions[j])]
            self.ggFiles += ['%s.gg_%s.txt' % (self.globalPrefix, self.controlRegions[j])]
            self.dggFiles += ['%s.dgg_%s.txt' % (self.globalPrefix, self.controlRegions[j])]
        for k in range(self.const.Nsplit):
            self.innerProductFiles += ['%s-%01d.inner_product.txt' % (self.globalPrefix, k)]
            self.ggFiles += ['%s-%01d.gg.txt' % (self.globalPrefix, k)]
            self.dggFiles += ['%s-%01d.dgg.txt' % (self.globalPrefix, k)]

        self.globalGradFiles, self.globalNormFiles, self.globalConjugateGradientFiles, self.globalControlSpaceFiles = [], [], [], []
        self.globalSensitivityFiles = []
        for k in range(self.const.NcontrolRegion):
            self.globalGradFiles += [self.globalPrefix + self.gradientFiles[k]]
            self.globalNormFiles += [self.globalPrefix + self.normFiles[k]]
            self.globalConjugateGradientFiles += [self.globalPrefix + self.conjugateGradientFiles[k]]
            self.globalControlSpaceFiles += [self.globalPrefix + self.controlForcingFiles[k]]
            self.globalSensitivityFiles += ['%s.adjoint_run_%s.txt' % (self.globalPrefix, self.controlRegions[k])]

        self.globalGradFiles += self.icGradientFiles
        self.globalControlSpaceFiles += self.icFiles
        for k in range(self.const.Nsplit):
            self.globalConjugateGradientFiles += ['%s-%01d.ic.conjugate_gradient.q' % (self.globalPrefix, k)]

        # previous gradient file
        self.previousGradFiles = [self.PREVDIR + '/previous.' + file for file in self.globalGradFiles]
        self.previousCGFiles = [self.PREVDIR + '/previous.' + file for file in self.globalConjugateGradientFiles]

        # Sort the files into subdirectories.
        self.innerProductFiles = [self.TXTDIR + '/' + file for file in self.innerProductFiles]
        self.ggFiles = [self.TXTDIR + '/' + file for file in self.ggFiles]
        self.dggFiles = [self.TXTDIR + '/' + file for file in self.dggFiles]

        self.diffFiles = [self.DIFFDIR + '/' + file for file in self.diffFiles]
        self.diffOutputFiles = [self.TXTDIR + '/' + file for file in self.diffOutputFiles]
        self.lagrangianFiles = [self.LGRNGDIR + '/' + file for file in self.lagrangianFiles]
        self.lagrangianOutputFiles = [self.TXTDIR + '/' + file for file in self.lagrangianOutputFiles]

        self.globalGradFiles = [self.GRADDIR + '/' + file for file in self.globalGradFiles]
        self.globalConjugateGradientFiles = [self.CGDIR + '/' + file for file in self.globalConjugateGradientFiles]
        self.globalSensitivityFiles = [self.TXTDIR + '/' + file for file in self.globalSensitivityFiles]

        self.baselineFiles = [self.BASELINEDIR + '/' + file for file in self.baselineFiles]

        return
