#!/bin/bash
#MOAB -l walltime=00:00:30:00,nodes=1:ppn=8
#MOAB -j oe
#MOAB -N pluto
#MOAB -q bbtest
#MOAB -m abe

cd "$PBS_O_WORKDIR"
module load apps/gcc/v4.8.4
module load apps/openmpi
mpirun -np 8 ./pluto -i pluto01.ini
