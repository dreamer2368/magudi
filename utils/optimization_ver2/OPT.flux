#!/bin/bash
#MSUB -l nodes=20
#MSUB -l partition=quartz
#MSUB -l walltime=3:00:00
#MSUB -m be
#MSUB -N AcousticMonopoleOptimization2
#MSUB -V
#MSUB -j oe -o result-%j.log
#MSUB -q pbatch

ppn=36
numProcs=$(($SLURM_NNODES*$ppn))
###numProcs=1056

export commandFile='AcousticMonopole.command.sh'
export decisionMakerCommandFile='AcousticMonopole.command.python.sh'
export nextDecisionMakerCommandFile='AcousticMonopole.command.python.ready.sh'

chmod u+x $commandFile

for k in {1..50}
do
    srun -N $SLURM_NNODES --mpi=none --mpibind=off flux start ./$commandFile
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

msub OPT.flux -l depend=$SLURM_JOBID
scontrol show job $SLURM_JOBID
