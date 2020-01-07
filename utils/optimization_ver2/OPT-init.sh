#!/bin/bash
#MSUB -l nodes=60
#MSUB -l partition=quartz
#MSUB -l walltime=0:30:00
#MSUB -m be
#MSUB -N SDML4
#MSUB -V
#MSUB -j oe -o result-%j.log
#MSUB -q pbatch

ppn=36
numProcs=$(($SLURM_NNODES*$ppn))
###numProcs=1056

export commandFile='WeiFreundSDML.command.sh'
export decisionMakerCommandFile='WeiFreundSDML.command.python.sh'
export nextDecisionMakerCommandFile='WeiFreundSDML.command.python.ready.sh'

bash initial.sh
if [ $? -ne 0 ]; then
    echo "initial step is not run successfully."
    exit -1
fi

#msub OPT.sh -l depend=$SLURM_JOBID
scontrol show job $SLURM_JOBID
