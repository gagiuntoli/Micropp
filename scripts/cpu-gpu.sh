#!/bin/bash

#SBATCH --job-name=test_parallel
#SBATCH -D .
#SBATCH --output=mpi_%j.out
#SBATCH --error=mpi_%j.err
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --gres=gpu:1

CPU_EXEC="../build-cpu/test/test3d_2"
GPU_EXEC="../build-cpu/test/test3d_2"

sizes=( 10 20 30 40 50 60 70 80 90 100 )

echo "SIZE    ASS_A     ASS_R      SOL      TOTAL" > times-cpu.txt

for i in ${sizes[@]}; do

	NN=$(( $i + 1 ))
	echo "CPU $i"
	fileo="cpu-${i}.txt"
	time $CPU_EXEC $NN 1 1 &> $fileo

	tassA=$(awk '/assembly_mat/{print $4}' $fileo)
	tassr=$(awk '/assembly_rhs/{print $4}' $fileo)
	tsol=$(awk '/ell_solve_cgpd/{print $4}' $fileo)

	echo "${i}x${i}x${i} $tassA $tassr $tsol $((tassA + tassr + tsol))" >> times-cpu.txt


	#echo "GPU $i"
	#fileo="gpu-${i}.txt"
	#time $GPU_EXEC $NN 1 1 &> $fileo

done
