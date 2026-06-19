#!/bin/bash
#MSUB -l nodes=160
#MSUB -l partition=quartz
#MSUB -l walltime=7:00:00
#MSUB -m be
#MSUB -N Kolmogorov
#MSUB -V
#MSUB -j oe -o result-%j.log
#MSUB -q pbatch

ppn=36
numProcs=$(($SLURM_NNODES*$ppn))
###numProcs=1056

export commandFile='KolmogorovFlow.command.sh'
export decisionMakerCommandFile='KolmogorovFlow.command.python.sh'
export nextDecisionMakerCommandFile='KolmogorovFlow.command.python.ready.sh'
export FLUX='/usr/global/tools/flux/toss_3_x86_64_ib/flux-0.17.0-pre-ft/bin/flux'
export LD_PRELOAD=/usr/global/tools/flux/toss_3_x86_64_ib/flux-0.17.0-pre-ft/lib/flux/libpmi2.so.0

chmod u+x $commandFile

for k in {1..10}
do
    mv $nextDecisionMakerCommandFile $decisionMakerCommandFile
    sh $decisionMakerCommandFile
    if [ $? -ne 0 ]; then
        echo "$decisionMakerCommandFile is not run successfully."
        exit -1
    fi
    srun -N $SLURM_NNODES --mpi=none --mpibind=off $FLUX start ./$commandFile
    if [ $? -ne 0 ]; then
        echo "$commandFile is not run successfully."
        exit -1
    fi
    scontrol show job $SLURM_JOBID
done

msub OPT.flux -l depend=$SLURM_JOBID
scontrol show job $SLURM_JOBID
