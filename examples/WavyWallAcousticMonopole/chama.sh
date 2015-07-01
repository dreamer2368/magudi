#!/bin/bash

NODES=$1
CORESPERNODE=8
TOTALCORES=$(($NODES * $CORESPERNODE))

export MAGUDI_MPIRUN="mpiexec --bind-to core --npernode ${CORESPERNODE} --n ${TOTALCORES}"

function setOption() {
    if grep -q "$1" magudi.inp
    then
	sed -i "s/^.*$1.*$/$1 = $2/g" magudi.inp
    else
	echo "$1 = $2" >> magudi.inp
    fi
}

function deleteOption() {
	if grep -q "$1" magudi.inp
    then
	sed -i "s/^.*$1.*$//g" magudi.inp
	fi
}

PREFIX=WavyWallAcousticMonopole
rm -f ${PREFIX}* Bootstrap/${PREFIX}*
python config.py

#get rid of the initial transient
cp ${PREFIX}.xyz magudi.inp ${PREFIX}_bc.dat Bootstrap/
cd Bootstrap
ln -sf /gscratch/dbuchta/magudi_target/magudi ./
NUMSTEPS=3600
setOption "number_of_timesteps" $NUMSTEPS
deleteOption "initial_condition_file"
deleteOption "control_mollifier_file"
deleteOption "target_mollifier_file"
deleteOption "mean_pressure_file"
setOption "disable_adjoint_solver" true
setOption "compute_time_average" true
${MAGUDI_MPIRUN} ../magudi

#computing the pressure average over the horizon of interest
cp ${PREFIX}-0000${NUMSTEPS}.q ../${PREFIX}.ic.q
cd ..
setOption "disable_adjoint_solver" true
setOption "compute_time_average" true
$MAGUDI_MPIRUN ./magudi

#using that average compute adjoint over that same horizon of interest
python -c "from config import *; mean_pressure(p3d.fromfile('${PREFIX}.mean.q')).save('${PREFIX}.mean_pressure.f')"
setOption "disable_adjoint_solver" false
setOption "compute_time_average" false
$MAGUDI_MPIRUN ./magudi
