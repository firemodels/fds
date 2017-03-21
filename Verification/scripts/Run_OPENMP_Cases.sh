#!/bin/bash
JOBPREFIX=O
OUTPUT=/tmp/openmp.out.$$
HOST=`hostname`
QUEUE=batch
FORCE=
KILL=
OOPT=
POPT=
benchbot_pid=~/.openmp_lock

#---------------------------------------------
#                   usage
#---------------------------------------------

function usage {
echo "Run_OPEN_Cases.sh [ -f -h -q queue_name -s"
echo "Runs OpenMP test cases"
echo ""
echo "Options"
echo "-f - force script to run"
echo "-h - display this message"
echo "-k - kill this script if it is running"
echo "-q queue_name - run cases using the queue queue_name"
echo "     default: batch"
echo "-s - stop FDS runs"
echo "-P - pass through -P option to qfds.sh"
echo "-O - pass through -O option to qfds.sh"
exit
}

#---------------------------------------------
#                   wait_cases_release_end
#---------------------------------------------

wait_cases_release_end()
{
   # Scans qstat and waits for cases to end
     while [[ `qstat -a | awk '{print $2 $4}' | grep $(whoami) | grep $JOBPREFIX` != '' ]]; do
        JOBS_REMAINING=`qstat -a | awk '{print $2 $4}' | grep $(whoami) | grep $JOBPREFIX | wc -l`
        echo "Waiting for ${JOBS_REMAINING} openmp benchmark cases to complete."
        sleep 30
     done
}

#---------------------------------------------
#                   LIST_DESCENDANTS
#---------------------------------------------

LIST_DESCENDANTS ()
{
  local children=$(ps -o pid= --ppid "$1")

  for pid in $children
  do
    LIST_DESCENDANTS "$pid"
  done

  echo "$children"
}

#---------------------------------------------
#                   main script
#---------------------------------------------

while getopts 'fhkO:P:q:s' OPTION
do
case $OPTION in
  f)
   FORCE=1
   ;;
  h)
   usage;
   ;;
  k)
   KILL=1;
   ;;
  q)
   QUEUE="$OPTARG"
   ;;
  O)
   OOPT="$OPTARG"
   ;;
  P)
   POPT="$OPTARG"
   ;;
  s)
   ./Run_FDS_Cases.sh -b -j $JOBPREFIX -s
   exit
   ;;
esac
done

if [ "$OOPT" != ]; then
  OOPT="-O $OOPT"
fi
if [ "$POPT" != ]; then
  POPT="-O $POPT"
fi

if [ "$KILL" == "1" ]; then
  pidrunning=0
  if [ -e $benchbot_pid ]; then
    PID=`head -1 $benchbot_pid`
    pidrunning=`ps -el |  awk '{print $4}'  | grep $PID | wc -l`
    rm -f $benchbot_pid
  fi
  if [ $pidrunning -gt 0 ]; then
    echo killing processes invoked by this script
    kill -9 $(LIST_DESCENDANTS $PID)
    echo "killing (PID=$PID)"
    kill -9 $PID
    JOBIDS=`qstat -a | grep $JOBPREFIX | awk -v user="$USER" '{if($2==user){print $1}}'`
    if [ "$JOBIDS" != "" ]; then
      echo killing jobs with Id: $JOBIDS
      qdel $JOBIDS
    fi
    echo process $PID killed
  else
    echo nothing to kill, Run_OPENMP_Cases.sh not running
  fi
  exit
fi

if [ "$FORCE" == "1" ]; then
  rm -f $benchbot_pid
fi
if [ -e $benchbot_pid ]; then
  echo "*** warning: another instance of this script is already running."
  echo "    If this is not the case, re-run this script with the -f option."
  exit 1
fi
echo $$ > $benchbot_pid

# run cases

./Run_FDS_Cases.sh -b -j $JOBPREFIX $OOPT $POPT -q $QUEUE
wait_cases_release_end

# copy 64x64x64 results to a file

echo $HOST
echo 64x64x64>$OUTPUT
for case in a b c d e f g h
do
CSV=../Timing_Benchmarks/openmp_test64${case}_devc.csv
time=`tail -1 $CSV | awk '{print $2}'`
echo $time>>$OUTPUT
done
echo "">>$OUTPUT
echo 128x128x128>>$OUTPUT

# copy 128x128x128 results to a file

for case in a b c d e f g h
do
CSV=../Timing_Benchmarks/openmp_test128${case}_devc.csv
time=`tail -1 $CSV | awk '{print $2}'`
echo $time>>$OUTPUT

# output or email results
done
if [ "x$EMAIL" == "x" ]; then
  cat $OUTPUT
else
  cat $OUTPUT | mail -s "openmp benchmark results on $HOST" $EMAIL
fi
rm $OUTPUT
rm $benchbot_pid
