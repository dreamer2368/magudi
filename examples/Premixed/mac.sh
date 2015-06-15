#!/bin/bash

export MAGUDI_MPIRUN="mpiexec -n 4"

function setOption() {
    if grep -q "$1" magudi.inp
    then
	sed -i "s/^.*$1.*$/$1 = $2/g" magudi.inp
    else
	echo "$1 = $2" >> magudi.inp
    fi
}

rm -f Premixed* Bootstrap/Premixed*
../../build/utils/premixed
cp Premixed.* Bootstrap
cd Bootstrap
$MAGUDI_MPIRUN ../../../build/magudi
cp Premixed-00001000.q ../Premixed.ic.q
cp Premixed-00001000.f ../Premixed.ic.f
cd ..
setOption "disable_adjoint_solver" true
setOption "compute_time_average" true
python config.py
$MAGUDI_MPIRUN ../../build/magudi
python -c "from config import *; mean_pressure(p3d.fromfile('Premixed.mean.q')).save('Premixed.mean_pressure.f')"
setOption "disable_adjoint_solver" false
setOption "compute_time_average" false
$MAGUDI_MPIRUN ../../build/magudi
