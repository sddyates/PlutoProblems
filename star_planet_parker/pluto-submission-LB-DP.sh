#!/bin/bash
#
#SBATCH --job-name=LB_DP
#SBATCH --output=LB_DP.txt
#
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --ntasks=32
#SBATCH --time=10-00:00:00
#SBATCH --mem-per-cpu=2000

module load openmpi
mpirun -np 32 ./pluto -i pluto_LB_DP.ini
