#!/bin/bash

if [ -d $2 ]; then
    FILENAME=MPI_openMP.log
else
    FILENAME=$3
fi
#rm -f $FILENAME

for i in {1..16}
do
    echo "Running on $i nodes"
    esubmit -n$i -t5 $1 $2 $3
    echo "Waiting for 15min more"
    sleep 15m # because hostfile is not job specific
done
