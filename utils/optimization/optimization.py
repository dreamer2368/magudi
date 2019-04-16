from base import *

NumProcs = 1
NumSearch = 10
prefix = ''
lineMinLog = prefix+'.line_minimization.txt'
CGLog = prefix+'.conjugate_gradient.txt'
dggFilename = prefix+'.dgg.txt'

forwardFilename = prefix+'.forward_run.txt'
adjointFilename = prefix+'.adjoint_run.txt'
controlForcingFilenames = [prefix+'.control_forcing.controlRegion.dat']
gradientFilenames = [prefix+'.gradient.controlRegion.dat']
CGFilenames = [prefix+'.conjugate_gradient.controlRegion.dat']
normFilenames = [prefix+'.norm.controlRegion.dat']

from mnbrak import *
from linmin import *
from frprmn import *

### conjugate-gradient procedure ###
###<------------------------------------------------------------- CURRENT STAGE: MOVE THIS WHENEVER GOING TO NEXT STAGE
### FORWARD RUN ###<==========================================================DON'T FORGET!!
### ADJOINT RUN ###<==========================================================DON'T FORGET!!

beforeLinmin(forwardFilename, adjointFilename,
             gradientFilenames, CGFilenames,
             normFilenames, initial=True)

### WHILE (REDUCTION>TOL):              ###


###     MNBRAK PROCEDURE                ###
setupInitialSteps(forwardFilename, CGFilenames,
                  controlForcingFilenames, zeroBaseline=True)
###     WHILE (NOT BRACKETED):          ###
###         NUMSEARCH FORWARD RUNS      ### <=================================DON'T FORGET!!
NextMnbrak(forwardFilename, CGFilenames,
           controlForcingFilenames, zeroBaseline=True)
###     MNBRAK PROCEDURE END            ###


###     LINMIN PROCEDURE                ###
###     WHILE (BRACKETSIZE > TOL):      ###
nextLinmin(forwardFilename, CGFilenames,
controlForcingFilenames, zeroBaseline=True, initial=True)
###         NUMSEARCH FORWARD RUNS      ### <=================================DON'T FORGET!!
###     LINMIN PROCEDURE END            ###

afterLinmin(forwardFilename, adjointFilename,
             gradientFilenames, CGFilenames,
             normFilenames, controlForcingFilenames)

###     GO BACK TO THE BEGINNING, OR CONJUGATE_GRADIENT PROCEDURE END       ###
