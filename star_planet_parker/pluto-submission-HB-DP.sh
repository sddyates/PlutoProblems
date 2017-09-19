#!/bin/bash
#
#SBATCH --job-name=HB_DP
#SBATCH --output=HB_DP.txt
#
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --ntasks=32
#SBATCH --time=10-00:00:00
#SBATCH --mem-per-cpu=2000

module load openmpi
mpirun -np 32 ./pluto -i pluto_HB_DP.ini
