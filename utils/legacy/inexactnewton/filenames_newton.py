from filenames import *         # inherit filenames from optimization_ver2

##################           optimization files              #####################

decisionMaker = 'inexactNewton.py'

# from os import path, getcwd
# ROOTDIR = getcwd()
# OUTDIR = '%s/out'%ROOTDIR
# DIFFDIR = '%s/diff'%ROOTDIR
# GRADDIR = '%s/grad'%ROOTDIR
# LGRNGDIR = '%s/lagrangian'%ROOTDIR
# TXTDIR = '%s/txt'%ROOTDIR
# CGDIR = '%s/cg'%ROOTDIR
# PREVDIR = '%s/previous'%ROOTDIR

# globalInputFile = 'magudi.inp'
# globalPrefix = 'AcousticMonopole'
# globalCommandFile = globalPrefix + '.command.sh'
# decisionMakerCommandFile = globalPrefix+'.command.python.ready.sh'
# lineMinLog = globalPrefix+'.line_minimization.txt'
# CGLog = globalPrefix+'.conjugate_gradient.txt'
# forwardLog = globalPrefix + '.cost_reduction.txt'
# forwardLogColumns
# gradientLog = globalPrefix + '.gg_log.txt'
# gradientLogColumns

RESDIR = '%s/residual'%ROOTDIR
NEWTONDIR = '%s/newton'%ROOTDIR
NEWTONADJDIR = '%s/newton-adjoint'%ROOTDIR
dirList += [RESDIR,NEWTONDIR,NEWTONADJDIR]

residualLog = globalPrefix + '.residual_log.txt'
residualLogColumns = ['total']
for k in range(Nsplit):
    residualLogColumns += ['residual-%d'%k]

##################           global level files              #####################

# controlRegions = ['controlRegion']
# gradientFiles, normFiles, conjugateGradientFiles, controlForcingFiles = [], [], [], []

##################           slice level files              #####################

# prefixes, inputFiles, outputFiles, commandFiles = [], [], [], []
# icFiles, diffFiles, diffOutputFiles             = [], [], []
# matchingForwardFiles, matchingAdjointFiles      = [], []
# icAdjointFiles, icGradientFiles, directories    = [], [], []
# lagrangianFiles, lagrangianOutputFiles          = [], []

##################           global+slice level files              #####################

# innerProductFiles, ggFiles, dggFiles = [], [], []
# globalGradFiles, globalNormFiles,
# globalConjugateGradientFiles, globalControlSpaceFiles = [], [], [], []
# globalSensitivityFiles = []

# previousGradFiles = [PREVDIR + '/previous.' + file for file in globalGradFiles]
# previousCGFiles = [PREVDIR + '/previous.' + file for file in globalConjugateGradientFiles]

# residual: residual of constraints between segments
# newton adjoint: adjoint newton step files, == number of constraints
# newton: newton step files, = number of contraints + number of control forcing files
newtonFiles = []
for j in range(NcontrolRegion):
    newtonFiles += ['%s.newton_%s.dat'%(globalPrefix,controlRegions[j])]

residualFiles, newtonAdjointFiles           = [], []
icLinearizedFiles, matchingLinearizedFiles  = [], []
diffResidualFiles                           = []
for k in range(Nsplit):
    residualFiles += ['%s-%01d.residual.q'%(globalPrefix,k)]
    newtonFiles += ['%s-%01d.newton.q'%(globalPrefix,k)]
    newtonAdjointFiles += ['%s-%01d.newton_adjoint.q'%(globalPrefix,k)]
    icLinearizedFiles += ['%s-%01d.ic.linearized.q'%(globalPrefix,k)]
    diffResidualFiles += ['%s-%01d.diff.residual.q'%(globalPrefix,k)]

    kOffset = startTimestep + Nts * (k+1)
    matchingLinearizedFiles += ['%s-%01d-%08d.linearized.q'%(globalPrefix,k,kOffset)]

newtonFiles = [ NEWTONDIR + '/' + file for file in newtonFiles ]
newtonAdjointFiles = [ NEWTONADJDIR + '/' + file for file in newtonAdjointFiles ]
residualFiles = [ RESDIR + '/' + file for file in residualFiles ]
diffResidualFiles = [ DIFFDIR + '/' + file for file in diffResidualFiles ]
