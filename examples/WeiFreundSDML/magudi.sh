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

function toggleOption() {
    if grep -q "$1 = false" magudi.inp
    then
	sed -i "s/$1 = false/$1 = true/g" magudi.inp
    elif grep -q "$1 = true" magudi.inp
    then
	sed -i "s/$1 = true/$1 = false/g" magudi.inp
    fi
}

rm -f WeiFreundSDML* Bootstrap/WeiFreundSDML*
python generateMesh.py
python generateTargetState.py
python generateInitialCondition.py
cp WeiFreundSDML.xyz WeiFreundSDML.target.q WeiFreundSDML.ic.q Bootstrap
cd Bootstrap
$MAGUDI_MPIRUN ../magudi
cp WeiFreundSDML-00033600.q ../WeiFreundSDML.ic.q
cd ..
setOption "disable_adjoint_solver" true
setOption "compute_time_average" true
$MAGUDI_MPIRUN ./magudi
python generateMeanPressure.py
python generateMollifiers.py
toggleOption "disable_adjoint_solver"
toggleOption "compute_time_average"
$MAGUDI_MPIRUN ./magudi
