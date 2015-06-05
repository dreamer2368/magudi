#!/bin/bash

export MAGUDI_MPIRUN="mpiexec --bind-to-core --npernode 8 --n 256"

function setOption() {
    if grep -q "$1" magudi.inp
    then
	sed -i "s/^.*$1.*$/$1 = $2/g" magudi.inp
    else
	echo "$1 = $2" >> magudi.inp
    fi
}

rm -f WeiFreundSDML* Bootstrap/WeiFreundSDML*
python config.py
cp WeiFreundSDML.xyz WeiFreundSDML.target.q WeiFreundSDML.ic.q Bootstrap
cd Bootstrap
$MAGUDI_MPIRUN ../magudi
cp WeiFreundSDML-00033600.q ../WeiFreundSDML.ic.q
cd ..
setOption "disable_adjoint_solver" true
setOption "compute_time_average" true
$MAGUDI_MPIRUN ./magudi
python -c "from config import *; mean_pressure(p3d.fromfile('WeiFreundSDML.mean.q')).save('WeiFreundSDML.mean_pressure.f')"
python generateMeanPressure.py
setOption "disable_adjoint_solver" false
setOption "compute_time_average" false
$MAGUDI_MPIRUN ./magudi
