import argparse

__all__ = ['parser']

descriptionString = 'Optimization script generator.\n\n'                                                          \
                    'Each time this produces a command shell file that will go through following loop:\n'       \
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
parser.add_argument('input_file', metavar='string', type=str,
                    help='filename for the optimization input file.\n')
parser.add_argument('--mode', metavar='string', action='store', dest='mode',
                    help='mode for the optimizer:\n'
                         '  setup: set up initial directories and files.\n'
                         '  schedule: write a command file for the next action.\n'
                         '  log_result: log the result for the command file.\n')
parser.add_argument('--result', metavar='integer', type=int,
                    help='log the result of the command file.\n')
