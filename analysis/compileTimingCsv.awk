#!/bin/awk -f
# Massage the .log files into a csv
# Example call:
# ./compileTimingCsv.awk  Naive_openMP.log Improved_openMP.log R-side_parallelism.log > test.csv
BEGIN{printf("Setting, cores, time\n")}
{$1 = FILENAME; gsub(".log","")}
    {$2 = FILENAME; gsub("_"," ")} #; printf "%s \t", $1}
    /cores/{printf "%s %s, %s, ", $1, $2, substr($0,length($0),1)}
    /elapsed/{getline; print $(NF)}
    
