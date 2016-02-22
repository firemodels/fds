#!/bin/bash

dir=$1
fdsfile=$2
curdir=`pwd`

SVNROOT=`pwd`/../..
cd $SVNROOT
SVNROOT=`pwd`
dir=$SVNROOT/Verification/$dir
csvoutfile=$SVNROOT/Utilities/Scripts/fds_timing_stats.csv
cd $curdir

fdsbase=${fdsfile%.*}
outfile=$fdsbase.out
cpufile=${fdsbase}_cpu.csv
cd $dir

#echo 'FDS Case,Wall Clock Time (s),CPU Time (s),Number of Cells,Number of Time Steps,Performance Metric (1e-6)' > $SVNROOT/Utilities/Scripts/fds_timing_stats.csv

   # Grep for wall clock time
   WALL_CLOCK_TIME_VALUE=`grep -H "Total Elapsed Wall Clock Time (s):" "$outfile" | awk -F' ' '{print $(NF)}'`

   # Grep for CPU time and units
   TOTAL_CPU_TIME=0.0
   CPU_TIME=`cat "$cpufile" | awk '{if(NR>1)print}' | awk -F',' '{print $(2)}'`
   for j in $CPU_TIME
   do
      jafter=`echo ${j} | sed -e 's/[eE]+*/\\*10\\^/'`
      jafter=`echo "$jafter" | bc`
      TOTAL_CPU_TIME=`echo "$TOTAL_CPU_TIME+$jafter" | bc`
   done

   # Grep for number of cells in each dimension
   X_CELLS=`grep -H "Cells in the X" $outfile | awk -F' ' '{print $NF}'`
   Y_CELLS=`grep -H "Cells in the Y" $outfile | awk -F' ' '{print $NF}'`
   Z_CELLS=`grep -H "Cells in the Z" $outfile | awk -F' ' '{print $NF}'`

   # Sum over the number of cells (for multi-mesh cases)
   for j in $X_CELLS
   do
      let NUM_TOTAL_CELLS+=`echo "$X_CELLS * $Y_CELLS * $Z_CELLS" | bc`
   done

   # Grep for number of time steps
   NUM_TIME_STEPS=`grep -H "Time Step     " $outfile | tail -n 1 | awk -F' ' '{print $(NF-4)}'`

   # Calculate nondimensional performance metric
   # Skip over cases with no time steps
   if [[ $NUM_TIME_STEPS -eq 0 ]]
   then
      NUM_TIME_STEPS=0
      PERFORMANCE=0
   else
      let PERFORMANCE=`echo "1000000*$TOTAL_CPU_TIME/($NUM_TOTAL_CELLS*$NUM_TIME_STEPS)" | bc`
   fi

   # Write results to fds_timing_stats.csv file
#   echo "$FILE,$WALL_CLOCK_TIME_VALUE,$CPU_TIME,$NUM_TOTAL_CELLS,$NUM_TIME_STEPS,$PERFORMANCE" >> $csvoutfile
echo WALL_CLOCK_TIME_VALUE=$WALL_CLOCK_TIME_VALUE
echo TOTAL_CPU_TIME=$TOTAL_CPU_TIME
echo X_CELLS=$X_CELLS
echo Y_CELLS=$Y_CELLS
echo Z_CELLS=$Z_CELLS
echo NUM_TIME_STEPS=$NUM_TIME_STEPS
echo PERFORMANCE=$PERFORMANCE

   cd $curdir
