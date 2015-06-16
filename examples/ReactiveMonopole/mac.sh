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

rm -f ReactiveMonopole*
../../build/utils/premixed
python generate_mollifiers.py
$MAGUDI_MPIRUN ../../build/magudi
