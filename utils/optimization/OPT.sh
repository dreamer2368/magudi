#!/bin/bash
#MSUB -l nodes=10
#MSUB -l partition=quartz
#MSUB -l walltime=1:00:00:00
#MSUB -m be
#MSUB -N AcousticMonopoleOptimization
#MSUB -V
#MSUB -j oe -o result-%j.log
#MSUB -q pbatch

ppn=36
numProcs=$(($SLURM_NNODES*$ppn))
###numProcs=1056

export commandFile='AcousticMonopole.command.sh'
export decisionMakerCommandFile='AcousticMonopole.command.python.sh'

flag=true
while [ $flag=true ]
do
    sh $commandFile
    sh $decisionMakerCommandFile
    scontrol show job $SLURM_JOBID
done

scontrol show job $SLURM_JOBID
