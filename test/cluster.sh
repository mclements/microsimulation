#!/bin/bash

if [ -d $2 ]; then
    FILENAME=RunSim_simple.log
else
    FILENAME=$2
fi
#rm -f $FILENAME

cd $HOME
. ~/bin/setEnv.sh
. ~/bin/createHostFile.sh
#mpirun -n 1 --hostfile hostfile R --slave -f \
#~/microsimulation/test/cluster_mic.R
#export OMP_NUM_THREADS=1
mpirun -n 1 --hostfile hostfile R --slave -f $1 >> $2
