#!/bin/sh
#PBS -N cuzr_bond_order
#PBS -l select=1:ncpus=1
#PBS -j oe
#PBS -o $PBS_O_WORKDIR
module load /public/software/mpi/openmpi/1.8.3/intel 
export OMP_NUM_THREADS=1
cd $PBS_O_WORKDIR
#mpirun -np 16 
lmp_mpi -in  in.orientorder
date
