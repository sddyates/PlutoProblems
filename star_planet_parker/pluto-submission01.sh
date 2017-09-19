#!/bin/bash
#MOAB -l walltime=10:00:00:00,nodes=2:ppn=16
#MOAB -j oe
#MOAB -N pluto
#MOAB -q bblong
#MOAB -m abe

cd "$PBS_O_WORKDIR"
module load apps/gcc/v4.8.4
module load apps/openmpi
mpirun -np 32 ./pluto -i pluto01.ini
