from parser import parser
from optimizer import Optimizer
from inputs import InputParser

if __name__ == '__main__':
    args = parser.parse_args()
    if ((args.mode == "log_result") and (args.result is None)):
        raise RuntimeError("Specify the result as '--result integer'!\n")

    config = InputParser(args.input_file)
    optim = Optimizer(config)

    if (args.mode == "schedule"):
        optim.schedule()
    elif (args.mode == "log_result"):
        optim.checkResult(args.result)
    elif (args.mode == "setup"):
        optim.setupDirectories()

    else:
        print("choose the action. for possible choices, add '-h' flag and run it again.")
