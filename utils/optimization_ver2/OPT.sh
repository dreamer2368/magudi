#!/bin/bash
#MSUB -l nodes=120
#MSUB -l partition=quartz
#MSUB -l walltime=1:10:00
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

###bash initial.sh
###if [ $? -ne 0 ]; then
###    echo "initial step is not run successfully."
###    exit -1
###fi

for k in {1..10}
do
    bash $commandFile
    if [ $? -ne 0 ]; then
        echo "$commandFile is not run successfully."
        exit -1
    fi
    mv $nextDecisionMakerCommandFile $decisionMakerCommandFile
    sh $decisionMakerCommandFile
    if [ $? -ne 0 ]; then
        echo "$decisionMakerCommandFile is not run successfully."
        exit -1
    fi
    scontrol show job $SLURM_JOBID
done

msub OPT.sh -l depend=$SLURM_JOBID
scontrol show job $SLURM_JOBID
