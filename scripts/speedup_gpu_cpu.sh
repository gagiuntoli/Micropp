#!/bin/bash

#SBATCH --job-name=test_parallel
#SBATCH -D .
#SBATCH --output=mpi_%j.out
#SBATCH --error=mpi_%j.err
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --gres=gpu:1

CPU_EXEC="../build_p9_cpu/test/test3d_1"
GPU_EXEC="../build_p9_gpu/test/test3d_1"

sizes=( 25 50 75 100 )

for i in ${sizes[@]}; do

	echo "CPU $i"
	time $CPU_EXEC $i 0 1 &> cpu-${i}.txt

	echo "GPU $i"
	time $GPU_EXEC $i 0 1 &> gpu-${i}.txt

done
