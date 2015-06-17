#!/bin/bash

export MAGUDI_MPIRUN="mpiexec --bind-to core --npernode 8 --n 16"

function setOption() {
    if grep -q "$1" magudi.inp
    then
	sed -i "s/^.*$1.*$/$1 = $2/g" magudi.inp
    else
	echo "$1 = $2" >> magudi.inp
    fi
}

rm -f AcousticMonopole* Bootstrap/AcousticMonopole*
python config.py
cp AcousticMonopole.xyz Bootstrap
cd Bootstrap
$MAGUDI_MPIRUN ../magudi
cp AcousticMonopole-00000240.q ../AcousticMonopole.ic.q
cd ..
setOption "disable_adjoint_solver" true
setOption "compute_time_average" true
$MAGUDI_MPIRUN ./magudi
python -c "from config import *; mean_pressure(p3d.fromfile('AcousticMonopole.mean.q')).save('AcousticMonopole.mean_pressure.f')"
setOption "disable_adjoint_solver" false
setOption "compute_time_average" false
$MAGUDI_MPIRUN ./magudi
