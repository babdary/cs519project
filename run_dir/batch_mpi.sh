#!/bin/bash
#

#SBATCH -p defq   # kuyruk ismi
##SBATCH -A abasol

#SBATCH --job-name=mpi_c

#SBATCH --time=100
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
##SBATCH --mem-per-cpu=80000

##SBATCH --mem=240000
cd /home/bbilgin/cs519project/fw-par/row_distributed
ulimit -s 10240
mpirun ./main inputMatrix.txt
cat $HOSTFILE
rm -f $HOSTFILE
