#!/bin/bash
cd $HOME
. ~/bin/setEnv.sh
. ~/bin/createHostFile.sh
#mpirun -n 1 --hostfile hostfile R --slave -f \
#~/src/ki/microsimulation/test/cluster_mic.R
#export OMP_NUM_THREADS=1
mpirun -n 1 --hostfile hostfile R --slave -f \
~/microsimulation/test/cluster_mic2.R
