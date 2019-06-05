#!/bin/bash

##mpic++ multi-gpu-mpi.cpp -o multi-gpu-mpi -L ../build-gpu/ -lmicropp -I ../include/ -acc -std=c++11

#SBATCH --job-name="micropp"
#SBATCH --workdir=.
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --exclusive
#SBATCH --ntasks=320
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:4
#SBATCH --qos=debug
#SBATCH --time=01:00:00

module purge
module load gcc/6.4.0 bsc/commands pgi/18.10 cuda/10.1 cmake/3.13.2 spectrum/10.2

declare -a hosts=(`scontrol show hostnames $SLURM_NODELIST`)

echo "${hosts[0]}" > ./my_hostfile
echo "${hosts[0]}" >> ./my_hostfile
echo "${hosts[0]}" >> ./my_hostfile
echo "${hosts[0]}" >> ./my_hostfile
echo "${hosts[1]}" >> ./my_hostfile
echo "${hosts[1]}" >> ./my_hostfile
echo "${hosts[1]}" >> ./my_hostfile
echo "${hosts[1]}" >> ./my_hostfile


N=80
NGP=256
STEPS=5

#EXEC="../build-gpu/test/test3d_2"
EXEC="../../test/multi-gpu-mpi"

export OMP_NUM_THREADS=1
export CUDA_VISIBLE_DEVICES="0,1,2,3"

#time mpirun -np 1 $EXEC 100 1 20
#time nvprof -o nvprof-profile-${N}-${NGP}.out $EXEC $N $NGP 20
#time $EXEC 50 4 20

rm -rf times.txt

for i in 8; do

	echo "running" $i MPI processes
	fileo="output-${N}-${NGP}-${i}.out"
	time mpirun -np $i --machinefile ./my_hostfile -mca rmaps seq $EXEC $N $NGP $STEPS > $fileo
	#time mpirun -np $i --machinefile ./my_hostfile -mca rmaps seq $EXEC $N $NGP $STEPS > $fileo
	tim=$(awk '/time =/{print $3}' $fileo)
	echo $i $tim >> times.txt

#	mpirun -np 8 --machinefile ./my_hostfile -mca rmaps seq ./your_MPI_binary

done
