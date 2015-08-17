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



DIR=("1em16" "1em8" "1em4" "1em3" "1em2" "1em1")
AMP=("1.e-16" "1.e-8" "1.e-4" "1.e-3" "1.e-2" "1.e-1")

#DIR=("1em2" "1em1")
#AMP=("1.e-2" "1.e-1")
for i in "${!DIR[@]}"; do 
     rm -rf ${DIR[$i]}
	mkdir ${DIR[$i]}
	cp ${PREFIX}* *.inp ${DIR[$i]}/
	cd ${DIR[$i]}
	ln -sf /gscratch/dbuchta/magudi_target/magudi ./
	setOption "cost_functional_wallPenalty" "${AMP[$i]}"
	${MAGUDI_MPIRUN} ../magudi
	cd ../
done


