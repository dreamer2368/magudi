from constants import *

##################           optimization files              #####################

decisionMaker = 'optimization.py'

globalPrefix = 'AcousticMonopole'

globalCommandFile = globalPrefix + '.command.sh'
decisionMakerCommandFile = globalPrefix+'.command.python.ready.sh'

lineMinLog = globalPrefix+'.line_minimization.txt'
CGLog = globalPrefix+'.conjugate_gradient.txt'

forwardLog = globalPrefix + '.cost_reduction.txt'
forwardLogColumns = ['total']
for k in range(Nsplit):
    forwardLogColumns += ['cost_functional%d'%k]
for k in range(Nsplit-1):
    forwardLogColumns += ['matching_penalty%d'%k]

gradientLog = globalPrefix + '.gg_log.txt'
gradientLogColumns = ['total']
for j in range(NcontrolRegion):
    gradientLogColumns += ['control_forcing%d'%j]
for k in range(Nsplit-1):
    gradientLogColumns += ['initial_condition%d'%(k+1)]

##################           global level files              #####################

controlRegions = ['controlRegion']
gradientFiles, normFiles, conjugateGradientFiles, controlForcingFiles = [], [], [], []
for k in range(NcontrolRegion):
    gradientFiles += ['.gradient_%s.dat'%controlRegions[k]]
    normFiles += ['.norm_%s.dat'%controlRegions[k]]
    conjugateGradientFiles += ['.conjugate_gradient_%s.dat'%controlRegions[k]]
    controlForcingFiles += ['.control_forcing_%s.dat'%controlRegions[k]]

##################           slice level files              #####################

prefixes, inputFiles, outputFiles, commandFiles = [], [], [], []
icFiles, diffFiles, diffOutputFiles             = [], [], []
matchingForwardFiles, matchingAdjointFiles      = [], []
icAdjointFiles, icGradientFiles, directories    = [], [], []
lagrangianFiles, lagrangianOutputFiles          = [], []
for k in range(Nsplit):
    prefixes += ['%s-%01d'%(globalPrefix,k)]
    inputFiles += ['magudi-%01d.inp'%(k)]
    outputFiles += ['%s-%01d.forward_run.txt'%(globalPrefix,k)]
    commandFiles += ['%s-%01d.command.sh'%(globalPrefix,k)]
    icFiles += ['%s-%01d.ic.q'%(globalPrefix,k)]

    kOffset = startTimestep + NtimestepOffset * (k+1)
    diffFiles += ['%s-%01d.diff.q'%(globalPrefix,k)]                             #use only 0, ... , Nsplit-2
    diffOutputFiles += ['%s-%01d.diff.txt'%(globalPrefix,k)]                     #use only 0, ... , Nsplit-2
    matchingForwardFiles += ['%s-%01d-%08d.q'%(globalPrefix,k,kOffset)]          #use only 0, ... , Nsplit-2
    matchingAdjointFiles += ['%s-%01d-%08d.adjoint.q'%(globalPrefix,k,kOffset)]  #use only 0, ... , Nsplit-2
    lagrangianFiles += ['%s-%01d.lagrangian_multiplier.q'%(globalPrefix,k)]      #use only 0, ... , Nsplit-2
    lagrangianOutputFiles += ['%s-%01d.lagrangian.txt'%(globalPrefix,k)]         #use only 0, ... , Nsplit-2

    kOffset = startTimestep + NtimestepOffset * k
    icAdjointFiles += ['%s-%01d-%08d.adjoint.q'%(globalPrefix,k,kOffset)]        #use only 1, ... , Nsplit-1
    icGradientFiles += ['%s-%01d.ic.adjoint.q'%(globalPrefix,k)]                 #use only 1, ... , Nsplit-1
    directories += ['%01d'%k]

diffFiles[-1] = ''
diffOutputFiles[-1] = ''
matchingForwardFiles[-1] = ''
matchingAdjointFiles[-1] = ''
lagrangianFiles[-1] = ''
lagrangianOutputFiles[-1] = ''
icAdjointFiles[0] = ''
icGradientFiles[0] = ''

##################           global+slice level files              #####################

innerProductFiles, ggFiles, dggFiles = [], [], []
for j in range(NcontrolRegion):
    innerProductFiles += ['%s.inner_product_%s.txt'%(globalPrefix,controlRegions[j])]
    ggFiles += ['%s.gg_%s.txt'%(globalPrefix,controlRegions[j])]
    dggFiles += ['%s.dgg_%s.txt'%(globalPrefix,controlRegions[j])]
for k in range(1,Nsplit):
    innerProductFiles += ['%s-%01d.inner_product.txt'%(globalPrefix,k)]
    ggFiles += ['%s-%01d.gg.txt'%(globalPrefix,k)]
    dggFiles += ['%s-%01d.dgg.txt'%(globalPrefix,k)]

globalGradFiles, globalNormFiles, globalConjugateGradientFiles, globalControlSpaceFiles = [], [], [], []
globalSensitivityFiles = []
for k in range(NcontrolRegion):
    globalGradFiles += [globalPrefix + gradientFiles[k]]
    globalNormFiles += [globalPrefix + normFiles[k]]
    globalConjugateGradientFiles += [globalPrefix + conjugateGradientFiles[k]]
    globalControlSpaceFiles += [globalPrefix + controlForcingFiles[k]]
    globalSensitivityFiles += ['%s.adjoint_run_%s.txt'%(globalPrefix,controlRegions[k])]

globalGradFiles += icGradientFiles[1:]
globalControlSpaceFiles += icFiles[1:]
for k in range(1,Nsplit):
    globalConjugateGradientFiles += ['%s-%01d.ic.conjugate_gradient.q'%(globalPrefix,k)]


# dggFilename = prefix+'.dgg.txt'
# commandFilename = prefix+'.command.sh'
