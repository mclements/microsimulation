#!/bin/bash

# Add dependencies
module add easy
module load R/3.0.2

# Create a hostfile specefying the available nodes
. CreateHostfile.sh

# Go to home, where the hostfile is located and possible output will be saved
cd $HOME

# Inititating one node with mpirun and the others from within R
mpirun -n 1 --hostfile hostfile R --slave -f \
~/src/ki/microsimulation/test/cluster_mic.R
