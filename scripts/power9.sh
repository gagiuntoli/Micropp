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

time ../build-gpu/test/test_damage 100 20
