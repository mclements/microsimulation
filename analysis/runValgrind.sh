#!/bin/bash
. ~/bin/setEnv.sh

# Calls RunSim.R for 1 & 8 cores. The result is saved in a callgrind file.

for NCORES in 1 8; do
    echo "Number of cores:"$NCORES
    export OMP_NUM_THREADS=$NCORES
    R --vanilla -d "valgrind --tool=callgrind" < RunSim.R
done
