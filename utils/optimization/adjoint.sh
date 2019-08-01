#!/bin/bash
#MSUB -l nodes=30
#MSUB -l partition=quartz
#MSUB -l walltime=12:00:00
#MSUB -N 20000timestep_ADJOINT
#MSUB -m be
#MSUB -V
#MSUB -j oe -o result-%j.log
#MSUB -q pbatch

ppn=36
###numProcs=$(($SLURM_NNODES*$ppn))
numProcs=1056

export MAGUDI_MPIRUN="srun"
export EXEC='forward'

function setOption() {
    if grep -q "$1" magudi.inp
    then
	sed -i "s/^.*$1.*$/$1 = $2/g" magudi.inp
    else
	echo "$1 = $2" >> magudi.inp
    fi
}

export FORWARD='MultiblockJet.forward_run.txt'
export ADJOINT='MultiblockJet.adjoint_run.txt'

srun -N $SLURM_NNODES -n $numProcs ./adjoint $ADJOINT
scontrol show job $SLURM_JOBID

