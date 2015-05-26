#!/bin/awk -f
# Massage the .log files into a csv
# Example call:
# ./compileTimingCsv.awk  Naive_openMP.log Improved_openMP.log R-side_parallelism.log > test.csv
# or
# ./compileTimingCsv.awk O2.log O3.log Ofast_funroll.log > test.csv
BEGIN{printf("Setting, cores, time\n")}
{fname = FILENAME}
{sub(/.log/, "", fname)}
{sub(/_/, " ", fname)}
/cores/{printf "%s, %s, ", fname, substr($0,length($0),1)}
/elapsed/{getline; print $(NF)}
