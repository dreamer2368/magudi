#!/bin/bash
#MSUB -l nodes=32:ppn=8
#MSUB -l walltime=10:00:00
#MSUB -V
#MSUB -j oe
#MSUB -o log.o%j

export MAGUDI_MPIRUN="srun"

function setOption() {
    if grep -q "$1" magudi.inp
    then
	sed -i "s/^.*$1.*$/$1 = $2/g" magudi.inp
    else
	echo "$1 = $2" >> magudi.inp
    fi
}

rm -f Cylinder* Bootstrap/Cylinder*
python config.py
cp Cylinder.xyz Cylinder.target.q Bootstrap
cd Bootstrap
$MAGUDI_MPIRUN ../magudi
cp Cylinder-01000000.q ../Cylinder.ic.q
cd ..
setOption "disable_adjoint_solver" true
setOption "compute_time_average" true
$MAGUDI_MPIRUN ./magudi
python -c "from config import *; mean_pressure(p3d.fromfile('Cylinder.mean.q')).save('Cylinder.mean_pressure.f')"
setOption "disable_adjoint_solver" false
setOption "compute_time_average" false
$MAGUDI_MPIRUN ./magudi
