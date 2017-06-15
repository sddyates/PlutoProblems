#!/bin/bash
#MOAB -l walltime=10:00:00:00,nodes=10:ppn=16
#MOAB -j oe
#MOAB -N pluto
#MOAB -q bblong
#MOAB -m abe

cd "$PBS_O_WORKDIR"
module load apps/gcc/v4.8.4
mpirun -np 160 ./pluto -i pluto_spherical.ini -restart 8
