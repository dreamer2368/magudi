#!/bin/bash
#MSUB -l nodes=30
#MSUB -l partition=quartz
#MSUB -l walltime=5:00:00
#MSUB -N 20000timestep_FORWARD
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

msub adjoint.sh -l depend=$SLURM_JOBID

export FORWARD='MultiblockJet.forward_run.txt'
export ADJOINT='MultiblockJet.adjoint_run.txt'

cp 12/MultiblockJet.forward_run.txt ./ 
cp 12/MultiblockJet.control_forcing_controlRegion.E.dat ./ 
mv MultiblockJet.gradient_controlRegion.E.dat previous.MultiblockJet.gradient_controlRegion.E.dat
mv MultiblockJet.conjugate_gradient_controlRegion.E.dat previous.MultiblockJet.conjugate_gradient_controlRegion.E.dat
cp 12/MultiblockJet.control_forcing_controlRegion.W.dat ./ 
mv MultiblockJet.gradient_controlRegion.W.dat previous.MultiblockJet.gradient_controlRegion.W.dat
mv MultiblockJet.conjugate_gradient_controlRegion.W.dat previous.MultiblockJet.conjugate_gradient_controlRegion.W.dat
cp 12/MultiblockJet.control_forcing_controlRegion.N.dat ./ 
mv MultiblockJet.gradient_controlRegion.N.dat previous.MultiblockJet.gradient_controlRegion.N.dat
mv MultiblockJet.conjugate_gradient_controlRegion.N.dat previous.MultiblockJet.conjugate_gradient_controlRegion.N.dat
cp 12/MultiblockJet.control_forcing_controlRegion.S.dat ./ 
mv MultiblockJet.gradient_controlRegion.S.dat previous.MultiblockJet.gradient_controlRegion.S.dat
mv MultiblockJet.conjugate_gradient_controlRegion.S.dat previous.MultiblockJet.conjugate_gradient_controlRegion.S.dat

srun -N $SLURM_NNODES -n $numProcs ./forward $FORWARD
scontrol show job $SLURM_JOBID

