#!/bin/bash
#MOAB -l walltime=2:00:00:00,nodes=1:ppn=16
#MOAB -j oe
#MOAB -N pluto
#MOAB -q bbdefault
#MOAB -m abe

cd "$PBS_O_WORKDIR"
module load apps/gcc/v4.8.4
module load apps/openmpi
mpirun -np 16 ./pluto -i pluto_HUV.ini
