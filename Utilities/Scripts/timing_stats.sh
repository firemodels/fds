#!/bin/bash
# this script assumes it is being run in Verification
is_benchmark="no"
while getopts 'Ad:p:o:tT:' OPTION
do
case $OPTION  in
  A)
   is_benchmark="yes"
   ;;
  d)
   dir="$OPTARG"
   ;;
  p)
   dummy=1
   ;;
  o)
   dummy="$OPTARG"
   ;;
  t)
   dummy=1
   ;;
  T)
   dummy=1
   ;;
esac
done
shift $(($OPTIND-1))

# if we are saving benchmark times and this is not a benchmark case then exit
# otherwise save the times

if [ "$save_benchmark" == "yes" ]; then
  if [ "$is_benchmark" == "no" ]; then
     exit
  fi
fi

fdsfile=$1
curdir=`pwd`

fdsbase=${fdsfile%.*}
outfile=$fdsbase.out
cpufile=${fdsbase}_cpu.csv
cd $dir
if ! [ -e $outfile ]; then
  outfile=${fdsbase}_cat.out
  if ! [ -e $outfile ]; then
    exit
  fi
fi
if ! [ -e $cpufile ]; then
  cpufile=${fdsbase}_cat_cpu.csv
  if ! [ -e $cpufile ]; then
    exit
  fi
fi
   # Grep for wall clock time
   WALL_CLOCK_TIME_VALUE=`grep -H "Total Elapsed Wall Clock Time (s):" "$outfile" | awk -F' ' '{print $(NF)}'`

   # Grep for CPU time and units
   TOTAL_CPU_TIME=0.0
   CPU_TIME=`cat "$cpufile" | awk '{if(NR>1)print}' | awk -F',' '{print $(NF)}'`
   for j in $CPU_TIME
   do
      jafter=`echo ${j} | sed 's/\([0-9.]\+\)[eE][+]\?\([0-9-]\+\)/\1*10^\2/g'`
      jafter=`echo "$jafter" | bc -l `
      TOTAL_CPU_TIME=`echo "$TOTAL_CPU_TIME+$jafter" | bc `
   done
   TOTAL_CPU_TIME=$(echo "$TIMING_CPU_TIME" | sed 's/0*$//; s/\.$//')
   if [ "$TOTAL_CPU_TIME" == "" ]; then
     TOTAL_CPU_TIME=0.0
   fi
   CPU_TIME=$TOTAL_CPU_TIME

   # Grep for number of cells in each dimension
   X_CELLS="`grep -H "Cells in the X" $outfile | awk -F' ' '{print $NF}'`"
   set -- $X_CELLS
   let numx=$#
   Y_CELLS=`grep -H "Cells in the Y" $outfile | awk -F' ' '{print $NF}'`
   Z_CELLS=`grep -H "Cells in the Z" $outfile | awk -F' ' '{print $NF}'`
   # Sum over the number of cells (for multi-mesh cases)
   NUM_TOTAL_CELLS=0
   for i in $(seq 1 $numx);do
       XI=`echo $X_CELLS | awk -v var="$i" '{print $var}'`
       YI=`echo $Y_CELLS | awk -v var="$i" '{print $var}'`
       ZI=`echo $Z_CELLS | awk -v var="$i" '{print $var}'`
       sumxyz=`echo "$XI * $YI * $ZI" | bc`
       NUM_TOTAL_CELLS=`echo "$NUM_TOTAL_CELLS + $sumxyz" | bc`
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
   echo "$fdsfile,$WALL_CLOCK_TIME_VALUE,$CPU_TIME,$NUM_TOTAL_CELLS,$NUM_TIME_STEPS,$PERFORMANCE"

   cd $curdir
