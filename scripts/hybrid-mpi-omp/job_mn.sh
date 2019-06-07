#!/bin/bash

#SBATCH --job-name="script"
#SBATCH --workdir=.
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --ntasks=4
#SBATCH --exclusive
#SBATCH --time=00:10:00
#SBATCH --qos=debug

export OMP_NUM_THREADS=4
export EXTRAE_CONFIG_FILE=./extrae.xml
export LD_PRELOAD=${EXTRAE_HOME}/lib/libompitrace.so

#export EXTRAE_CONFIG_FILE=extrae.xml
#export LD_PRELOAD=${EXTRAE_HOME}/lib/libmpitrace.so

EXEC="../../test/hybrid-mpi-omp 20 12 10"

time mpirun -np 4 $EXEC
