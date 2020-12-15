#!/bin/bash
#MSUB -l nodes=160
#MSUB -l partition=quartz
#MSUB -l walltime=6:00:00
#MSUB -m be
#MSUB -N Kolmogorov-initial
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

chmod u+x initial-forward.sh
chmod u+x initial-adjoint.sh

srun -N $SLURM_NNODES --mpi=none --mpibind=off $FLUX start ./initial-forward.sh
if [ $? -ne 0 ]; then
    echo "initial forward is not run successfully."
    exit -1
fi
scontrol show job $SLURM_JOBID

srun -N $SLURM_NNODES --mpi=none --mpibind=off $FLUX start ./initial-adjoint.sh
if [ $? -ne 0 ]; then
    echo "initial adjoint is not run successfully."
    exit -1
fi
scontrol show job $SLURM_JOBID

###cat <<EOF > KolmogorovFlow.command.python.ready.sh
###python3 optimization.py 1 -initial_cg -zero_baseline
###EOF
cat <<EOF > KolmogorovFlow.command.python.ready.sh
python3 optimization.py 1 -initial_cg
EOF

msub OPT.flux -l depend=$SLURM_JOBID
scontrol show job $SLURM_JOBID
