#!/bin/bash

export JOBPREFIX=MP_
OUTFILE=openmp_summary.csv
QFDS="../../Utilities/Scripts/qfds.sh  -I  -t -q bench1 " 
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
  ./makecase64.sh 64 test64a$arg
  $QFDS -o 1 test64a$arg.fds

  ./makecase64.sh 64 test64d$arg
  $QFDS -o 4 test64d$arg.fds

  ./makecase128.sh 128 test128a$arg
  $QFDS -o 1 test128a$arg.fds

  ./makecase128.sh 128 test128d$arg
  $QFDS -o 4 test128d$arg.fds
done

wait_cases_end

echo 64 1 thread, Host, 64 4 threads, Host, 128 1 thread, Host, 128 4 threads, Host>$OUTFILE
for i in `seq 1 $NCASES`; do
arg=0$i
  if [ $i -gt 9 ]; then
    arg=$i
  fi
  
  HOST1=`grep Host test64a$arg.log  | awk '{print $2}'`
  HOST2=`grep Host test64d$arg.log  | awk '{print $2}'`
  HOST3=`grep Host test128a$arg.log | awk '{print $2}'`
  HOST4=`grep Host test128d$arg.log | awk '{print $2}'`
  TIME1=`grep Time test64a$arg.out  | grep Stepping | awk '{print $7}'`
  TIME2=`grep Time test64d$arg.out  | grep Stepping | awk '{print $7}'`
  TIME3=`grep Time test128a$arg.out | grep Stepping | awk '{print $7}'`
  TIME4=`grep Time test128d$arg.out | grep Stepping | awk '{print $7}'`
  echo "$TIME1,$HOST1,$TIME2,$HOST2,$TIME3,$HOST3,$TIME4,$HOST4">>$OUTFILE
done
echo complete
