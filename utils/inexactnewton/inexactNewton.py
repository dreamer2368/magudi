from base import *
from newton import *
import subprocess
import argparse

from os import path
for dir in dirList:
    if (not path.exists(dir)):
        subprocess.check_call('mkdir -p %s'%dir, shell=True)
if (not path.exists('linminLog')):
    subprocess.check_call('mkdir -p linminLog', shell=True)

descriptionString = 'Inexact Newton script generator.\n\n'                                                      \
                    'Each time this produces a command shell file that will go through following loop:\n'       \
                    ' forward\n'                                                                                \
                    ' residual\n'                                                                               \
                    ' while(cg optimization continues):\n'                                                      \
                    '     adjoint\n'                                                                            \
                    '     linearized\n'                                                                         \
                    '     linmin\n'                                                                             \
                    '     cgstep\n'                                                                             \
                    ' newtonstep\n'

parser = argparse.ArgumentParser(description = descriptionString,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('procedure', metavar='string', choices=['forward','residual','adjoint',
                                                            'linearized','cgstep','newtonstep'],
                    help='choose the optimization action. Choices are:\n'
                         ' forward\n'
                         ' residual\n'
                         ' adjoint\n'
                         ' linearized\n'
                         ' linmin\n'
                         ' cgstep\n'
                         ' newtonstep')
parser.add_argument('-zero_control', metavar='boolean',
                    action='store_const', const=True, default=False,
                    help='True when not using control paramters')
# parser.add_argument('-initial_cg', metavar='boolean',
#                     action='store_const', const=True, default=False,
#                     help='use at the first loop of cg optimization.')
# parser.add_argument('-linmin_initial', metavar='boolean',
#                     action='store_const', const=True, default=False,
#                     help='use at the first loop of line minimization.')
# parser.add_argument('-linmin_stop', metavar='boolean',
#                     action='store_const', const=True, default=False,
#                     help='use when force-stopping the line minimization before reaching the tolerance.')

args = parser.parse_args()
zeroControl = args.zero_control

switcher = {
    'forward': forward,
    'residual': residual,
    'adjoint': adjoint,
    'linearized': linearized,
    'linmin': linmin,
    'cgstep': cgstep,
    'newtonstep': newtonstep
}

action = switcher.get(args.procedure, lambda: 'Invalid procedure')
success = action(zeroControl)
if (not success):
    print ('%s procedure is failed.'%args.procedure)
