#!/bin/bash
#SBATCH --nodes=60
# The partititon option is not supported, please log onto the
# Machine you want to use before submitting your job.
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dreamer2368@gmail.com
#SBATCH --job-name=MultiblockJet-check
#SBATCH --output=result-%j.log
#SBATCH --partition=pbatch

ppn=36
numProcs=$(($SLURM_NNODES*$ppn))
###numProcs=1056

export FLUX='/usr/global/tools/flux/toss_3_x86_64_ib/default/bin/flux'
###export FLUX='/usr/global/tools/flux/toss_3_x86_64_ib/flux-0.17.0-pre-ft/bin/flux'
###export LD_PRELOAD=/usr/global/tools/flux/toss_3_x86_64_ib/flux-0.17.0-pre-ft/lib/flux/libpmi2.so.0
export J0file='MultiblockJet.forward.0.txt'
export gradfile='MultiblockJet.adjoint.txt'
export commandFile='MultiblockJet.command.sh'

index=0

for k in {1..1}
do
    python3 checkGradientAccuracy.py optim.yml --mode schedule
    if [ $? -ne 0 ]; then
        echo "Scheduling is not run successfully."
        exit -1
    fi
    chmod u+x ${commandFile}
    srun -N $SLURM_NNODES --mpi=none --mpibind=off $FLUX start ./${commandFile}
    export RESULT=$?
    if [ $RESULT -eq 1 ]; then
      echo "Gradient check finished."
      exit 0
    else
      echo $RESULT
      python3 checkGradientAccuracy.py optim.yml --mode log_result --result $RESULT
    fi
done

###python3 -c "from checkGradientAccuracy import writeStepForwardRunCommand; writeStepForwardRunCommand(filename='step-forward.sh', idx=$index);"
###if [ $? -ne 0 ]; then
###    echo "step forward script is not created successfully."
###    exit -1
###fi
###
###chmod u+x step-forward.sh
###
###srun -N $SLURM_NNODES --mpi=none --mpibind=off $FLUX start ./step-forward.sh
###if [ $? -ne 0 ]; then
###    echo "step forward is not run successfully."
###    exit -1
###fi
###scontrol show job $SLURM_JOBID
###
###python3 -c "from checkGradientAccuracy import computeGradientAccuracy; computeGradientAccuracy(filename='MultiblockJet.gradient_accuracy.txt', idx=$index, '$J0file', '$gradfile');"
###if [ $? -ne 0 ]; then
###    echo "computeGradientAccuracy is not run successfully."
###    exit -1
###fi
###scontrol show job $SLURM_JOBID

sbatch check.flux -l depend=$SLURM_JOBID
scontrol show job $SLURM_JOBID
##MSUB #-V # This is default behavior in SLURM.
##MSUB #-j # Removing -e line (if it exists) , slurm will combine.
