#!/bin/bash

export JOBPREFIX=MP_
OUTFILE=openmp_summary.csv
QFDS="../../Utilities/Scripts/qfds.sh  -I  -t -q firebot " 
NCASES=5

#---------------------------------------------
#                   wait_cases_end
#---------------------------------------------

wait_cases_end()
{
   # Scans qstat and waits for cases to end
     while [[ `qstat -a | awk '{print $2 $4}' | grep $(whoami) | grep $JOBPREFIX` != '' ]]; do
        JOBS_REMAINING=`qstat -a | awk '{print $2 $4}' | grep $(whoami) | grep $JOBPREFIX | wc -l`
        echo "Waiting for ${JOBS_REMAINING} openmp benchmark cases to complete."
        sleep 30
     done
}

for i in `seq 1 $NCASES`; do
arg=0$i
  if [ $i -gt 9 ]; then
    arg=$i
  fi
  ./makecase64.sh 64 openmp_test64a$arg
  $QFDS -o 1 openmp_test64a$arg.fds

  ./makecase64.sh 64 openmp_test64d$arg
  $QFDS -o 4 openmp_test64d$arg.fds

  ./makecase128.sh 128 openmp_test128a$arg
  $QFDS -o 1 openmp_test128a$arg.fds

  ./makecase128.sh 128 openmp_test128d$arg
  $QFDS -o 4 openmp_test128d$arg.fds
done

wait_cases_end

echo 64 1 thread, 64 4 threads, 128 1 thread, 128 4 threads>$OUTFILE
for i in `seq 1 $NCASES`; do
arg=0$i
  if [ $i -gt 9 ]; then
    arg=$i
  fi
  TIME1=`grep Time openmp_test64a$arg.out | grep Stepping | awk '{print $7}'`
  TIME2=`grep Time openmp_test64d$arg.out | grep Stepping | awk '{print $7}'`
  TIME3=`grep Time openmp_test128a$arg.out | grep Stepping | awk '{print $7}'`
  TIME4=`grep Time openmp_test128d$arg.out | grep Stepping | awk '{print $7}'`
  echo "$TIME1,$TIME2,$TIME3,$TIME4">>$OUTFILE
done
echo complete
