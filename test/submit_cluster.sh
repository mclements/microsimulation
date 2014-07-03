#!/bin/bash
qsub -q fas_high -l nodes=4:ppn=8,walltime=00:02:00 cluster.sh
