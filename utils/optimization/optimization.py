from base import *
from mnbrak import *
from linmin import *
from frprmn import *
import subprocess
import argparse

descriptionString = 'Optimization decision maker.\n\n'                                                          \
                    'Each time this produces a command shell file that will go through each steps in the following loop:\n'       \
                    ' while(cg optimization continues):\n'                                                      \
                    '     forward run\n'                                                                        \
                    '     adjoint run\n'                                                                        \
                    '     beforeLinmin\n\n'                                                                     \
                    '     setupInitialSteps\n'                                                                  \
                    '     while(finding initial bracket):\n'                                                    \
                    '         intermediate forward runs\n'                                                      \
                    '         nextMnbrak\n'                                                                     \
                    '     while(narrowing bracket enough)\n'                                                    \
                    '         nextLinmin\n'                                                                     \
                    '         intermediate forward runs\n'                                                      \
                    '     afterLinmin'

parser = argparse.ArgumentParser(description = descriptionString,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('procedure', metavar='integer', type=int, nargs=1, default=[0],
                    help='choose the optimization action. Choices are:\n'
                         ' 1: beforeLinmin\n'
                         ' 2: mnbrak - setupInitialSteps\n'
                         ' 3: mnbrak - nextMnbrak\n'
                         ' 4: linmin - nextLinmin\n'
                         ' 5: afterLinmin')
parser.add_argument('-zero_baseline', metavar='boolean',
                    action='store_const', const=True, default=False,
                    help='use when current control forcing is zero')
parser.add_argument('-initial_cg', metavar='boolean',
                    action='store_const', const=True, default=False,
                    help='use at the first loop of cg optimization.')
parser.add_argument('-linmin_initial', metavar='boolean',
                    action='store_const', const=True, default=False,
                    help='use at the first loop of line minimization.')
parser.add_argument('-linmin_stop', metavar='boolean',
                    action='store_const', const=True, default=False,
                    help='use when force-stopping the line minimization before reaching the tolerance.')

args = parser.parse_args()
zeroBaseline = args.zero_baseline
initial = args.initial_cg
linminInitial = args.linmin_initial
linminStop = args.linmin_stop

action = args.procedure[0]

if( action==1 ):
    success = beforeLinmin(forwardFilename, adjointFilename,
                           gradientFilenames, CGFilenames,
                           normFilenames, initial, zeroBaseline)
    if (success!=0):
        print ('pre-processing for line minimization is failed.')

elif( action==2 ):
    success = setupInitialSteps(forwardFilename, CGFilenames,
                                controlForcingFilenames, zeroBaseline)
    if (success!=0):
        print ('initialization for mnbrak is failed.')

elif( action==3 ):
    retry = False
    case = NextMnbrak(forwardFilename, CGFilenames,
                      controlForcingFilenames, zeroBaseline, retry)
    if (case==0):
        bracketed = True
    elif (case>0):
        bracketed = False
    else:
        print ('mnbrak is failed.')

elif( action==4 ):
    case = nextLinmin(forwardFilename, CGFilenames,
                      controlForcingFilenames, zeroBaseline, initial=linminInitial, stop=linminStop)
    if (case==0):
        linminContinue = False

elif( action==5 ):
    case = afterLinmin(forwardFilename, adjointFilename,
                       gradientFilenames, CGFilenames,
                       normFilenames, controlForcingFilenames)
    if (case==0):
        CG_continue = False

else:
    print("choose the action. for possible choices, add '-h' flag and run it again.")
