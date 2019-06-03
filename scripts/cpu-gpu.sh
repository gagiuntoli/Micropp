#!/bin/bash

#SBATCH --job-name=test_parallel
#SBATCH -D .
#SBATCH --output=mpi_%j.out
#SBATCH --error=mpi_%j.err
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --gres=gpu:1

CPU_EXEC="../build-cpu/test/test3d_1"
GPU_EXEC="../build-cpu/test/test3d_1"

sizes=( 25 50 75 100 )

for i in ${sizes[@]}; do

	echo "CPU $i"
	fileo="cpu-${i}.txt"
	time $CPU_EXEC $i 0 1 &> $fileo

	echo "GPU $i"
	fileo="gpu-${i}.txt"
	time $GPU_EXEC $i 0 1 &> $fileo

done
