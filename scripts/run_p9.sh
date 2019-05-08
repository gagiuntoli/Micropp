#!/bin/bash

#SBATCH --job-name="test"
#SBATCH --workdir=.
#SBATCH --output=test_%j.out
#SBATCH --error=test_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=160
### #SBATCH --nodes=1
### #SBATCH --ntasks-per-node=4
### #SBATCH --qos=bsc_case
####SBATCH --qos=debug
#SBATCH --gres=gpu:4
#SBATCH --exclusive
#SBATCH --time=00:05:00

EXEC=""

export OMP_NUM_THREADS=1
export CUDA_VISIBLE_DEVICES="0,1,2,3"

time mpirun -np 1 $EXEC
