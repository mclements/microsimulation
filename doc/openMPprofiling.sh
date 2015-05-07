#!/bin/bash Rscript

# Calls RunSim.R for 1,2,4 & 8 cores. The result is saved in a log file where "elapsed time" is the relevant measure.

if [ -d $1 ]; then
    FILENAME=RunSim.log
else
    FILENAME=$1
fi
rm -f $FILENAME
for NCORES in 1 2 4 8; do
    echo "Number of cores:"$NCORES
    export OMP_NUM_THREADS=$NCORES
    echo "Number of cores:"$NCORES >> $FILENAME
    Rscript RunSim.R >> $FILENAME
done
