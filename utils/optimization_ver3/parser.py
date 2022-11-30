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
