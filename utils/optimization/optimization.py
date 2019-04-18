from base import *
from mnbrak import *
from linmin import *
from frprmn import *
import subprocess

zeroBaseline = False
initial = True

### conjugate-gradient procedure ###
CG_continue = True
while (CG_continue):
###<------------------------------------------------------------- CURRENT STAGE: MOVE THIS WHENEVER GOING TO NEXT STAGE
### FORWARD RUN ###<==========================================================DON'T FORGET!!
    subprocess.check_call('srun -n '+str(NumProcs)+' ./forward',shell=True)
### ADJOINT RUN ###<==========================================================DON'T FORGET!!
    subprocess.check_call('srun -n '+str(NumProcs)+' ./adjoint',shell=True)

    success = beforeLinmin(forwardFilename, adjointFilename,
                            gradientFilenames, CGFilenames,
                            normFilenames, initial)
    if (success!=0):
        print ('pre-processing for line minimization is failed.')
        break
    initial = False

### WHILE (REDUCTION>TOL):              ###


###     MNBRAK PROCEDURE                ###
    success = setupInitialSteps(forwardFilename, CGFilenames,
                                controlForcingFilenames, zeroBaseline)
    if (success!=0):
        print ('initialization for mnbrak is failed.')
        break
###     WHILE (NOT BRACKETED):          ###
    bracketed = False
    while (not bracketed):
###         NUMSEARCH FORWARD RUNS      ### <=================================DON'T FORGET!!
        for k in range(1,NumSearch+1):
            subprocess.check_call('cd '+str(k),shell=True)
            subprocess.check_call('srun -n '+str(NumProcs)+' ./forward',shell=True)
            subprocess.check_call('cd ..',shell=True)
            
#subprocess.check_call('sh intermediate_forward_runs.sh',shell=True)

        case = NextMnbrak(forwardFilename, CGFilenames,
                          controlForcingFilenames, zeroBaseline)
        if (case==0):
            bracketed = True
        elif (case>0):
            bracketed = False
        else:
            print ('mnbrak is failed.')
            break
            break           
###     MNBRAK PROCEDURE END            ###


###     LINMIN PROCEDURE                ###
    linminContinue = True
    linminInitial = False
    while (linminContinue):
###     WHILE (BRACKETSIZE > TOL):      ###
        case = nextLinmin(forwardFilename, CGFilenames,
                          controlForcingFilenames, zeroBaseline, initial=linminInitial)
        linminInitial = False
        if (case==0):
            linminContinue = False
            break
###         NUMSEARCH FORWARD RUNS      ### <=================================DON'T FORGET!!
        for k in range(1,NumSearch+1):
            subprocess.check_call('cd '+str(k),shell=True)
            subprocess.check_call('srun -n '+str(NumProcs)+' ./forward',shell=True)
            subprocess.check_call('cd ..',shell=True)
#subprocess.check_call('sh intermediate_forward_runs.sh',shell=True)
###     LINMIN PROCEDURE END            ###

    case = afterLinmin(forwardFilename, adjointFilename,
                       gradientFilenames, CGFilenames,
                       normFilenames, controlForcingFilenames)
    if (case==0):
        CG_continue = False

###     GO BACK TO THE BEGINNING, OR CONJUGATE_GRADIENT PROCEDURE END       ###
