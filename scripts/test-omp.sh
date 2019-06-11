#!/bin/bash

#SBATCH --job-name=OPENMP
#SBATCH --workdir=.
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --exclusive
#SBATCH --time=02:00:00
##SBATCH --qos=bsc_case
#SBATCH --qos=debug

N=30
NGP=480
STEPS=2

procs=(65 70 80 85 90 95 100)

rm -f times-${N}-${NGP}-${STEPS}.txt

for i in ${procs[@]}; do

	echo "Executing with $i Threads"
	export OMP_NUM_THREADS=${i}

	folder="res_${N}_${NGP}"
	mkdir -p ${folder}

	time ../build/test/test_omp $N $NGP $STEPS > ${folder}/res_${i}.dat
	tim=$(awk '/time =/{print $3}' ${folder}/res_${i}.dat)
	echo "$i $tim" >> times-${N}-${NGP}-${STEPS}.txt

done
