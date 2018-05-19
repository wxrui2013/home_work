#!/bin/sh
#PBS -N cuzr_s1_2000
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -o $PBS_O_WORKDIR
module load /public/software/mpi/openmpi/1.8.3/intel 
export OMP_NUM_THREADS=1
cd $PBS_O_WORKDIR
#mpirun -np 8 
lmp_mpi -in  gr.in
date
