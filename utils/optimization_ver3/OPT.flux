#!/bin/bash
#SBATCH --nodes=160
# The partititon option is not supported, please log onto the
# Machine you want to use before submitting your job.
#SBATCH --time=7:00:00
#SBATCH --mail-type=BEGIN,END
#SBATCH --job-name=Kolmogorov
#SBATCH --output=result-%j.log
#SBATCH --partition=pbatch

ppn=36
numProcs=$(($SLURM_NNODES*$ppn))
###numProcs=1056

export commandFile='KolmogorovFlow.command.sh'
export FLUX='/usr/global/tools/flux/toss_3_x86_64_ib/flux-0.17.0-pre-ft/bin/flux'
export LD_PRELOAD=/usr/global/tools/flux/toss_3_x86_64_ib/flux-0.17.0-pre-ft/lib/flux/libpmi2.so.0

chmod u+x $commandFile

for k in {1..10}
do
  python3 optimization.py optim.yml --mode schedule
  if [ $? -ne 0 ]; then
      echo "Scheduling is not run successfully."
      exit -1
  fi

  srun -N $SLURM_NNODES --mpi=none --mpibind=off $FLUX start ./$commandFile
  export RESULT=$?
  python3 optimization.py optim.yml --mode log_result --result $RESULT

  scontrol show job $SLURM_JOBID
done

sbatch OPT.flux -l depend=$SLURM_JOBID
scontrol show job $SLURM_JOBID
##MSUB #-V # This is default behavior in SLURM.
##MSUB #-j # Removing -e line (if it exists) , slurm will combine.
