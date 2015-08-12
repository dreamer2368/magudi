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
rm -f ${PREFIX}* 
python config.py

DIR=1em8
AMP=1.e-8
mkdir $DIR
cp ${PREFIX}* *.inp $DIR/
cd $DIR
ln -sf /gscratch/dbuchta/magudi_target/magudi ./
setOption "cost_functional_wallPenalty" "${AMP}"
${MAGUDI_MPIRUN} ../magudi
cd ../

DIR=1em4
AMP=1.e-4
mkdir $DIR
cp ${PREFIX}* *.inp $DIR/
cd $DIR
ln -sf /gscratch/dbuchta/magudi_target/magudi ./
setOption "cost_functional_wallPenalty" "${AMP}"
${MAGUDI_MPIRUN} ../magudi
cd ../


DIR=1em2
AMP=1.e-2
mkdir $DIR
cp ${PREFIX}* *.inp $DIR/
cd $DIR
ln -sf /gscratch/dbuchta/magudi_target/magudi ./
setOption "cost_functional_wallPenalty" "${AMP}"
${MAGUDI_MPIRUN} ../magudi
cd ../

DIR=1em1
AMP=1.e-1
mkdir $DIR
cp ${PREFIX}* *.inp $DIR/
cd $DIR
ln -sf /gscratch/dbuchta/magudi_target/magudi ./
setOption "cost_functional_wallPenalty" "${AMP}"
${MAGUDI_MPIRUN} ../magudi
cd ../

