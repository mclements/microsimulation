#!/bin/bash
. ~/bin/setEnv.sh

# Calls RunSim.R for 1,2,4 & 8 cores. The result is saved in a log file where "elapsed time" is the relevant measure.

if [ -d $1 ]; then
    FILENAME=RunSim_simple.log
else
    FILENAME=$1
fi
rm -f $FILENAME
export OMP_NUM_THREADS=1
Rscript RunSim_simple.R >> $FILENAME
