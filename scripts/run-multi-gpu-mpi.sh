#!/bin/bash

#SBATCH --job-name="micropp"
#SBATCH --workdir=.
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:4
#SBATCH --exclusive
#SBATCH --time=01:30:00

N=5
NGP=5
STEPS=5
N_MPI=5

EXEC="../build/multi-gpu-mpi"

export OMP_NUM_THREADS=1
export CUDA_VISIBLE_DEVICES="0,1,2,3"

#time mpirun -np 1 $EXEC 100 1 20
#time nvprof -o nvprof-profile-${N}-${NGP}.out $EXEC $N $NGP 20
#time $EXEC 50 4 20

rm -rf times.txt

for (( i=1; i<=$N_MPI; i++ )); do

	echo "running" $i MPI processes
	mpirun -np $i $EXEC $N $NGP $STEPS > output-${N}-${NGP}-${i}.out
	tim=$(awk '/time =/{print $3}' output-${N}-${NGP}-${i}.out)
	echo $i $tim >> times.txt

done
