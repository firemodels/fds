#!/bin/bash

# fds_timing_stats.sh
# Kristopher Overholt
# 8/24/2012

# This script checks the .out files for the FDS Verification Suite
# and generates a .csv file of the CPU time called fds_timing_stats.csv

SVNROOT=`pwd`/../..
cd $SVNROOT/Verification

echo 'FDS Case, CPU Time' > $SVNROOT/Utilities/Scripts/fds_timing_stats.csv
for i in */*.out
do
   FILE=$i
   CPU_TIME=`grep -H "Total CPU:" "$i" | tail -n 1 | awk -F' ' '{print $(NF-1), $NF}'`
   echo "$FILE, $CPU_TIME" >> $SVNROOT/Utilities/Scripts/fds_timing_stats.csv
done