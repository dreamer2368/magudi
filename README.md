# `magudi` tutorial

## Prerequisite
1. Fortran compiler: `ifort` or `gfortran`
2. MPI compiler: `mvapich` or `openmpi`

MPI compiler must be installed with the Fortran compiler that compiles `magudi`.

## Installation
Say the source code directory (mostly `magudi`) is named as `<magudi-src>`. We recommend making a separate directory for the installation, which we denote `<magudi-build>`. The installation is then executed via the following commands on the terminal,
```
mkdir <magudi-build>
cd <magudi-build>
cmake -D<cmake-var>=% <magudi-src>
make
```

The flag `-D<cmake-var>=%` is optional, to manually set `<cmake-var>` to `%`. Some useful flags `<cmake-var>` are:
* `CMAKE_BUILD_TYPE`: set to either `release` or `debug`
* `CMAKE_Fortran_COMPILER`: set to a specific fortran compiler binary file

## Example: acoustic monopole
### Configuration
We copy the example directory `<magudi-src>/examples/AcousticMonopole` to any location you like.
```
cp -r <magudi-src>/examples/AcousticMonopole/* ./
```

`config.py` file includes all routines to generate grid, function, and flow solution files in PLOT3D format. This needs a python file from `<magudi-src>/utils/python`.
```
cp <magudi-src>/utils/python/plot3dnasa.py ./
```

Note that `plot3dnasa.py` and `config.py` are written in python2. They are currently not compatible with python3. You can either import `config.py` and use the routines, or run `config.py` itself:

```
python2 config.py
```

This will generate:

* Grid file: `AcousticMonopole.xyz`
* Initial condition file: `AcousticMonopole.ic.q`
* Function files: `AcousticMonopole.target_mollifier.f`, `AcousticMonopole.control_mollifier.f`

### Baseline simulation
After compiling, main binary executables are generated in `<magudi-build>/bin/`. There are also auxiliary executables in `<magudi-build>/utils/`. Check the codes to see what tasks they do.

For baseline simulation, make a link to `<magudi-build>/bin/forward`:
```
ln -s <magudi-build>/bin/forward ./
```

`forward` requires two input files:
* `magudi.inp`: this includes all flags/parameters for simulation.
* `bc.dat`: this includes boundary condition information on grids. you can specify your own boundary condition file in `magudi.inp`.

To execute the baseline simulation, run:
```
./forward --output J0.txt
```

This saves the quantity-of-interest (QoI) J in `J0.txt`. Additionally, solutions at designated timesteps are saved. These will be required for adjoint simulation. Please check the code to see more optional arguments for `forward`.

### Adjoint simulation
For baseline simulation, make a link to `<magudi-build>/bin/adjoint`:
```
ln -s <magudi-build>/bin/adjoint ./
```

To execute the adjoint simulation, run:
```
./adjoint --output Grad0.txt
```

This saves the gradient magnitude in `Grad0.txt`. Gradient vector (forcing) will be saved as `AcousticMonopole.gradient.controlRegion.dat`. Additionally, solutions at designated timesteps are saved. These will be required for adjoint simulation. Please check the code to see more optional arguments for `adjoint`.

### Gradient accuracy check by finite-difference
We check the gradient accuracy by applying the control forcing along the gradient direction. In `magudi.inp`, change the `control_forcing_switch` flag,
```
controller_switch = true
```

To make the control forcing file, make a link to `<magudi-build>/bin/zaxpy`:
```
ln -s <magudi-build>/bin/zaxpy ./
```

`zaxpy` refers to `z=a*x+y`, where `z`, `x`, and `y` are the vector `.dat` files with the same length of the gradient file, and `a` is a scalar. We make a control forcing that is 0.0001 times smaller than the gradient,
```
./zaxpy AcousticMonopole.control_forcing.controlRegion.dat 1e-4 AcousticMonopole.gradient.controlRegion.dat`
```
where the argument `y` is automatically taken to be zero. Run another forward run with this file,
```
./forward --output J1.txt
```
where QoI is saved in `J1.txt`. A python3 routine to compute the finite-difference and check the accuracy is
```
import numpy as np
fID = open('J0.txt','r')
J0 = float(fID.read())
fID.close()
fID = open('Grad0.txt','r')
Grad0 = float(fID.read())
fID.close()
fID = open('J1.txt','r')
J1 = float(fID.read())
fID.close()

A1 = 1.0E-4
Grad1 = (J1-J0) / A1
Error = abs( (Grad1 - Grad0)/Grad0 )
print ("{:.16E}".format(A1), "{:.16E}".format(Error))
```

Following python3 subroutine does this job for multiple amplitudes and save the errors in `AcousticMonopole.gradient_accuracy.txt`:
```
import numpy as np
import subprocess

fID = open('J0.txt','r')
QoI0 = float(fID.read())
fID.close()
fID = open('Grad0.txt','r')
Grad0 = float(fID.read())
fID.close()

Nk = 20
Ak = 10.0**(-2.0-0.25*np.array(range(Nk)))
QoIk = np.zeros((Nk,),dtype=np.double)
Gradk = np.zeros((Nk,),dtype=np.double)
ek = np.zeros((Nk,),dtype=np.double)

for k in range(Nk):
    amp = Ak[k]
    command = ''
    command += './zaxpy AcousticMonopole.control_forcing.controlRegion.dat %.16E AcousticMonopole.gradient.controlRegion.dat`' % amp
    command += './forward --output J1.txt'
    fID = open('test-step.sh','w')
    fID.write(command)
    fID.close()
    subprocess.check_call('bash test-step.sh', shell=True)

    fID = open('J1.txt','r')
    QoIk[k] = float(fID.read())
    fID.close()

    Gradk[k] = (QoIk[k]-QoI0)/Ak[k]
    ek[k] = abs( (Gradk[k]-Grad0)/Grad0 )
    print ("{:.16E}".format(Ak[k]), "{:.16E}".format(QoIk[k]), "{:.16E}".format(Gradk[k]), "{:.16E}".format(ek[k]))

    fId = open(globalPrefix+'.gradient_accuracy.txt','a+')
    fId.write('%.16E\t%.16E\t%.16E\t%.16E\n' % (Ak[k], QoIk[k], Gradk[k], ek[k]))
    fId.close()
```
