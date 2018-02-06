#!/bin/bash

curdir=`pwd`
SVNROOT=`pwd`/../..
cd $SVNROOT
SVNROOT=`pwd`

tmpout=/tmp/timings.$$

# Mode argument: smokebot, or anything else (including no argument) for normal mode
MODE=$1
TIMING=$2

echo 'FDS Case,Wall Clock Time (s),CPU Time (s),Number of Cells,Number of Time Steps,Performance Metric (1e-6)' > $tmpout

export QFDS=$SVNROOT/Utilities/Scripts/timing_stats.sh
export RUNCFAST=$SVNROOT/Utilities/Scripts/timing_stats.sh
export RUNTFDS=$SVNROOT/Utilities/Scripts/timing_stats.sh
if [[ "$TIMING" == "" ]]; then 
  export save_benchmark=
else
  export save_benchmark="yes"
fi
cd $SVNROOT/Verification
if [[ "$MODE" == "smokebot" ]]; then
  scripts/SMV_Cases.sh >> $tmpout
  scripts/GEOM_Cases.sh >> $tmpout
  scripts/WUI_Cases.sh >> $tmpout
else
  ./FDS_Cases.sh >> $tmpout
  ./FDS_Benchmark_Cases_single.sh >> $tmpout
fi

TOTAL_CPU_TIME=0.0
CPU_TIME=`cat "$tmpout" | awk '{if(NR>1)print}' | awk -F',' '{print $(3)}'`
for j in $CPU_TIME
do
   jafter=`echo ${j} | sed -e 's/[eE]+*/\\*10\\^/'`
   jafter=`echo "$jafter" | bc`
   TOTAL_CPU_TIME=`echo "$TOTAL_CPU_TIME+$jafter" | bc`
done
echo `git rev-parse --short HEAD` >> $tmpout
echo $TOTAL_CPU_TIME >> $tmpout

cat $tmpout
rm $tmpout
export only_benchmark_times=
