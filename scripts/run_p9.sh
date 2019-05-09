#!/bin/bash

#SBATCH --job-name="micropp"
#SBATCH --workdir=.
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=160
#SBATCH --gres=gpu:4
#SBATCH --exclusive
#SBATCH --time=00:05:00

### #SBATCH --nodes=1
### #SBATCH --ntasks-per-node=4
### #SBATCH --qos=bsc_case
### #SBATCH --qos=debug

N=50
NGP=1

EXEC="../build-gpu/test/test3d_2"

export OMP_NUM_THREADS=4
export CUDA_VISIBLE_DEVICES="0,1,2,3"

#time mpirun -np 1 $EXEC 100 1 20
time nvprof -o nvprof-profile-${N}-${NGP}.out $EXEC $N $NGP 20
#time $EXEC 50 4 20
