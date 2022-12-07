#!/bin/bash
#SBATCH --nodes=160
# The partititon option is not supported, please log onto the
# Machine you want to use before submitting your job.
#SBATCH --time=7:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dreamer2368@gmail.com
#SBATCH --job-name=Kolmogorov
#SBATCH --output=result-%j.log
#SBATCH --partition=pbatch

ppn=36
numProcs=$(($SLURM_NNODES*$ppn))
###numProcs=1056

export commandFile='WeiFreundSDML.command.sh'

for k in {1..10}
do
  python3 optimization.py optim.yml --mode schedule
  if [ $? -ne 0 ]; then
      echo "Scheduling is not run successfully."
      exit -1
  fi

  bash $commandFile
  export RESULT=$?
  python3 optimization.py optim.yml --mode log_result --result $RESULT

  scontrol show job $SLURM_JOBID
done

sbatch OPT.sh -l depend=$SLURM_JOBID
scontrol show job $SLURM_JOBID
