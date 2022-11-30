from parser import parser
from optimizer import Optimizer
from inputs import InputParser

if __name__ == '__main__':
    args = parser.parse_args()
    zeroBaseline = args.zero_baseline
    initial = args.initial_cg
    linminInitial = args.linmin_initial
    linminStop = args.linmin_stop

    action = args.procedure[0]

    config = InputParser(args.input_file)
    optim = Optimizer(config)

    if (action == 1):
        success = optim.beforeLinmin(initial, zeroBaseline)
        if (success!=0):
            print ('pre-processing for line minimization is failed.')

    elif (action == 2):
        success = optim.setupInitialSteps(zeroBaseline)
        if (success != 0):
            print ('initialization for mnbrak is failed.')

    elif (action == 3):
        case = optim.nextMnbrak(zeroBaseline)
        if (case == 0):
            bracketed = True
        elif (case > 0):
            bracketed = False
        else:
            print ('mnbrak is failed.')

    elif (action == 4):
        case = optim.nextLinmin(zeroBaseline, initial=linminInitial, stop=linminStop)
        if (case == 0):
            linminContinue = False

    elif (action == 5):
        case = optim.afterLinmin(zeroBaseline)
        if (case == 0):
            CG_continue = False

    elif (action == -1):
        from os import path
        import subprocess
        for dir in optim.fl.dirList:
            if (not path.exists(dir)):
                subprocess.check_call('mkdir -p %s'%dir, shell=True)

    else:
        print("choose the action. for possible choices, add '-h' flag and run it again.")
