#!/bin/bash

# fds_timing_stats.sh
# Kristopher Overholt
# 8/24/2012

# This script checks the .out files for the FDS Verification Suite
# and generates a .csv file of the CPU time called fds_timing_stats.csv

# cd to Verification directory
SVNROOT=`pwd`/../..
cd $SVNROOT/Verification

# Write header information to fds_timing_stats.csv file
echo 'FDS Case,CPU Time (s),Number of Cells,Number of Time Steps,Performance Metric' > $SVNROOT/Utilities/Scripts/fds_timing_stats.csv

# Loop over all .out files in the Verification directory
for i in */*.out
do
   # Get current file name
   FILE=$i

   # Grep for CPU time and units
   CPU_TIME_VALUE=`grep -H "Total CPU:" "$i" | tail -n 1 | awk -F' ' '{print $(NF-1)}'`
   CPU_TIME_UNITS=`grep -H "Total CPU:" "$i" | tail -n 1 | awk -F' ' '{print $NF}'`

   # Convert min to s
   if [[ $CPU_TIME_UNITS = 'min' ]]
   then
      CPU_TIME=`echo "$CPU_TIME_VALUE*60" | bc`
   fi

   # Convert hr to s
   if [[ $CPU_TIME_UNITS = 'h' ]]
   then
      CPU_TIME=`echo "$CPU_TIME_VALUE*3600" | bc`
   fi

   # Grep for number of cells in each dimension
   X_CELLS=`grep -H "Cells in the X" "$i" | awk -F' ' '{print $NF}'`
   Y_CELLS=`grep -H "Cells in the Y" "$i" | awk -F' ' '{print $NF}'`
   Z_CELLS=`grep -H "Cells in the Z" "$i" | awk -F' ' '{print $NF}'`

   # Sum over the number of cells (for multi-mesh cases)
   for j in $X_CELLS
   do
      let NUM_TOTAL_CELLS+=`echo "$X_CELLS * $Y_CELLS * $Z_CELLS" | bc`
   done

   # Grep for number of time steps
   NUM_TIME_STEPS=`grep -H "Time Step     " "$i" | tail -n 1 | awk -F' ' '{print $(NF-4)}'`

   # Calculate nondimensional performance metric
   # Skip over cases with no time steps
   if [[ $NUM_TIME_STEPS -eq 0 ]]
   then
      NUM_TIME_STEPS=0
      PERFORMANCE=0
   else
      let PERFORMANCE=`echo "($CPU_TIME*$NUM_TOTAL_CELLS)/$NUM_TIME_STEPS" | bc`
   fi

   # Write results to fds_timing_stats.csv file
   echo "$FILE,$CPU_TIME,$NUM_TOTAL_CELLS,$NUM_TIME_STEPS,$PERFORMANCE" >> $SVNROOT/Utilities/Scripts/fds_timing_stats.csv

   # Reset variables
   NUM_TOTAL_CELLS=0
done
