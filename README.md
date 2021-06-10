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

This saves the quantity-of-interest J in `J0.txt`. Additionally, solutions at designated timesteps are saved. These will be required for adjoint simulation. Please check the code to see more optional arguments for `forward`.

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

`zaxpy` refers to $$\mathbf{z}=a\mathbf{x}+y$$
