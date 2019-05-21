#!/bin/bash

#SBATCH --job-name=test_parallel
#SBATCH -D .
#SBATCH --output=mpi_%j.out
#SBATCH --error=mpi_%j.err
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --gres=gpu:1

GPU_EXEC="../build_p9_gpu/test/test3d_4"

for ntask in {1..8}; do

	echo "OMP TASKS = $ntask"
	export OMP_NUM_THREADS=$ntask
	time $GPU_EXEC 50 24 1 &> speedup-multi-gpus-${ntask}.txt

done
