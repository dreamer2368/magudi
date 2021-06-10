# `magudi` tutorial

## Prerequisite
1. Fortran compiler: `ifort` or `gfortran`
2. MPI compiler: `mvapich` or `openmpi`

MPI compiler must be installed with the Fortran compiler that compiles `magudi`.

## Installation
Say the source code directory (mostly `magudi`) is named as `<magudi-src>`. We recommend making a separate directory for the installation, which we denote `<magudi-build>`. The installation is then executed via the following commands on the terminal,
>`mkdir <magudi-build>`
>
>`cd <magudi-build>`
>
>`cmake -D<cmake-var>=% <magudi-src>`
>
>`make`

The flag `-D<cmake-var>=%` is optional, to manually set `<cmake-var>` to `%`. Some useful flags `<cmake-var>` are:
* `CMAKE_BUILD_TYPE`: set to either `release` or `debug`
* `CMAKE_Fortran_COMPILER`: set to a specific fortran compiler binary file

## Example: acoustic monopole
We copy the example directory `<magudi-src>/examples/AcousticMonopole` to any location you like.
>`cp -r <magudi-src>/examples/AcousticMonopole/* ./`

`config.py` file includes all routines to generate grid, function, and flow solution files in PLOT3D format. This needs a python file from `<magudi-src>/utils/python`.

>`cp <magudi-src>/utils/python/plot3dnasa.py ./`

Note that `plot3dnasa.py` and `config.py` are written in python2. They are currently not compatible with python3. You can either import `config.py` and use the routines, or run `config.py` itself:

>`python2 config.py`

This will generate:

* Grid file: `AcousticMonopole.xyz`
* Random solution files: `AcousticMonopole-0.ic.q`, `AcousticMonopole-1.ic.q`
* Function files: `AcousticMonopole.target_mollifier.f`, `AcousticMonopole.control_mollifier.f`
