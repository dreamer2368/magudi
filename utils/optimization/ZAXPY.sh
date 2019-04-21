#!/bin/bash
#MSUB -l nodes=10
#MSUB -l partition=quartz
#MSUB -l walltime=0:10:00
#MSUB -N ZAXPY
#MSUB -m e
#MSUB -V
#MSUB -j oe -o result-%j.log
#MSUB -q pbatch

ppn=36
numProcs=$(($SLURM_NNODES*$ppn))
###numProcs=1056

echo srun -n $numProcs ./zaxpy $1 -$2 $3 $4
srun -n $numProcs ./zaxpy $1 -$2 $3 $4


scontrol show job $SLURM_JOBID

