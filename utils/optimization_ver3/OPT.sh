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

msub OPT.sh -l depend=$SLURM_JOBID
scontrol show job $SLURM_JOBID
