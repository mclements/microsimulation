#!/bin/bash
#PBS -o CLUSTER
#PBS -j oe
#PBS -m abe
#PBS -M mark.clements@ki.se
module load Apps/R/3.0.2
module load Rpkgs/DOSNOW
module load Rpkgs/RMPI
module add Rpkgs/RCPP/1.11.1
cd $PBS_O_WORKDIR
mpirun -n 1 R --slave -f cluster_mic.R
