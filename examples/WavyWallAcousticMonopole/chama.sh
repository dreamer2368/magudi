#!/bin/bash

NODES=$1
CORESPERNODE=16
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
cp ${PREFIX}.xyz magudi.inp ${PREFIX}.ambient_pressure.f ${PREFIX}_bc.dat Bootstrap/
cd Bootstrap
ln -sf /gscratch/dbuchta/magudi_target/magudi ./
NUMSTEPS=500
setOption "number_of_timesteps" $NUMSTEPS
deleteOption "initial_condition_file"
deleteOption "control_mollifier_file"
deleteOption "target_mollifier_file"
setOption "disable_adjoint_solver" true
setOption "compute_time_average" false
${MAGUDI_MPIRUN} ../magudi

echo 'ran magudi'

#computing the pressure average over the horizon of interest
cp ${PREFIX}-00000${NUMSTEPS}.q ../${PREFIX}.ic.q
cd ..

#setOption "disable_adjoint_solver" true
#setOption "compute_time_average" true
#$MAGUDI_MPIRUN ./magudi

#using that average compute adjoint over that same horizon of interest
#instead of the mean let's write an ambient pressure file
#python -c "from config import *; mean_pressure(p3d.fromfile('${PREFIX}.mean.q')).save('${PREFIX}.mean_pressure.f')"
#python -c "from config import *; ambient_pressure(p3d.fromfile('${PREFIX}.ambient.q')).save('${PREFIX}.ambient_pressure.f')"

setOption "disable_adjoint_solver" false
setOption "compute_time_average" false
NUMSTEPS=100
setOption "number_of_timesteps" $NUMSTEPS
$MAGUDI_MPIRUN ./magudi
