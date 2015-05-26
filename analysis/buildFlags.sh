#!/bin/bash
. ~/bin/setEnv.sh

# Wrinting flags to ~/.R/Makevars and then calls RunSim.R for 1,2,4 & 8 cores. The result is saved in a log file where "elapsed time" is the relevant output.
FLAGS=("-O2" "-O3" "-Ofast -funroll-all-loops")
FNAME=("O2.log" "O3.log" "Ofast_funroll.log")
echo ${FLAGS[0]}
for i in 0 1 2; do
    echo "Writing CXXFLAGS = "${FLAGS[i]}" in ~/.R/Makevars before building."
    echo "CXXFLAGS = "${FLAGS[i]} > ~/.R/Makevars
    R CMD INSTALL ~/microsimulation --preclean
    echo "Writing profiling results to: "${FNAME[i]}
    rm -f  ${FNAME[i]}
    for NCORES in 1 2 4 8; do
    	echo "Number of cores:"$NCORES
    	export OMP_NUM_THREADS=$NCORES
    	echo "Number of cores:"$NCORES >> ${FNAME[i]}
    	Rscript RunSim.R >> ${FNAME[i]}
    done
done

echo "Clearing ~/.R/Makevars"
rm -f ~/.R/Makevars
