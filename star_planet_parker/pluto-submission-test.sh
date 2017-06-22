#!/bin/bash
#
#SBATCH --job-name=pluto_SPI_test
#SBATCH --output=pluto_SPI_test.txt
#
#SBATCH --get-user-env
#SBATCH --export=NONE
#SBATCH --ntasks=8
#SBATCH --time=00-00:30:00
#SBATCH --mem-per-cpu=2000

module load openmpi
mpirun -np 8 ./pluto -i pluto_HUV.ini
