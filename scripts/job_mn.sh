#!/bin/bash
#SBATCH --job-name="script"
#SBATCH --workdir=.
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --qos=debug

./porcentages.sh
