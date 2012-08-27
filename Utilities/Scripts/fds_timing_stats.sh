#!/bin/bash

# fds_timing_stats.sh
# Kristopher Overholt
# 8/24/2012

# This script checks the .out files for the FDS Verification Suite
# and generates a .csv file of the CPU time called fds_timing_stats.csv

SVNROOT=`pwd`/../..
cd $SVNROOT/Verification

echo 'FDS Case,CPU Time,Number of Cells,Number of Time Steps' > $SVNROOT/Utilities/Scripts/fds_timing_stats.csv
for i in */*.out
do
   FILE=$i
   CPU_TIME=`grep -H "Total CPU:" "$i" | tail -n 1 | awk -F' ' '{print $(NF-1), $NF}'`
   X_CELLS=`grep -H "Cells in the X" "$i" | awk -F' ' '{print $NF}'`
   Y_CELLS=`grep -H "Cells in the Y" "$i" | awk -F' ' '{print $NF}'`
   Z_CELLS=`grep -H "Cells in the Z" "$i" | awk -F' ' '{print $NF}'`
   NUM_TIME_STEPS=`grep -H "Time Step     " "$i" | tail -n 1 | awk -F' ' '{print $(NF-4)}'`
   NUM_TOTAL_CELLS=0
   for j in $X_CELLS
   do
       let NUM_TOTAL_CELLS+=`echo "$X_CELLS * $Y_CELLS * $Z_CELLS" | bc`
   done
   echo "$FILE,$CPU_TIME,$NUM_TOTAL_CELLS,$NUM_TIME_STEPS" >> $SVNROOT/Utilities/Scripts/fds_timing_stats.csv
done