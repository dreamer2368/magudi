#!/bin/bash
#MSUB -l nodes=1
#MSUB -l partition=quartz
#MSUB -l walltime=00:30:00
### MSUB -m be
#MSUB -V
#MSUB -j oe -o result-%j.log
#MSUB -q pdebug

srun -n 30 pvbatch adjoint.py
