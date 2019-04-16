#!/bin/bash

#SBATCH --job-name=test_k80
#SBATCH --workdir=.
#SBATCH --output=k80_%j.out
#SBATCH --error=k80_%j.err
#SBATCH --ntasks=1
#SBATCH --gres gpu:1
#SBATCH --cpus-per-task=1
#SBATCH --constraint=k80
#SBATCH --time=00:02:00

EXEC="../build_mt/test/test_ell_mvp_openacc"
EXEC_ACC="../build_mt_gpu/test/test_ell_mvp_openacc"
#EXEC_ACC="pgprof -f -o a.prof ../build_mt_gpu/test/test_ell_mvp_openacc"
#EXEC="../build_mt/test/test3d_1"
#EXEC_ACC="../build_mt_gpu/test/test3d_1"
ARGS="100 20"

#time srun $EXEC 100 100
echo "Exec Normal Case"
time $EXEC $ARGS

echo 
echo "Exec OpenACC Case"
time $EXEC_ACC $ARGS
