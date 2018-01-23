#!/bin/bash

export JOBPREFIX=MP_
OUTFILE=openmp_summary.csv
QFDS="../../Utilities/Scripts/qfds.sh  -P -I  -t -q batch " 
NCASES=99

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
  base=t64${arg}
  ./makecase.sh 64 ${base}a
  ./makecase.sh 64 ${base}d
  $QFDS $base

  base=t128${arg}
  ./makecase.sh 128 ${base}a
  ./makecase.sh 128 ${base}d
  $QFDS $base
done

wait_cases_end

echo 64 1 thread, Host, 64 4 threads, Host, 128 1 thread, Host, 128 4 threads, Host>$OUTFILE
for i in `seq 1 $NCASES`; do
arg=0$i
  if [ $i -gt 9 ]; then
    arg=$i
  fi
  arglog=${arg} 
  arga=${arg}a
  argd=${arg}d
  HOST1=
  HOST2=
  HOST3=
  HOST4=
  TIME1=
  TIME2=
  TIME3=
  TIME4=
  if [ -e t64$arglog.log ]; then
    HOST1=`grep Host t64$arglog.log  | tail -1 | awk '{print $2}'`
    HOST2=`grep Host t64$arglog.log  | tail -1 | awk '{print $2}'`
  fi
  if [ -e t128$arglog.log ]; then
    HOST3=`grep Host t128$arglog.log | tail -1 | awk '{print $2}'`
    HOST4=`grep Host t128$arglog.log | tail -1 | awk '{print $2}'`
  fi
  if [ -e t64$arga.out ]; then
    TIME1=`grep Time t64$arga.out  | grep Stepping | awk '{print $7}'`
  fi
  if [ -e t64$argd.out ]; then
    TIME2=`grep Time t64$argd.out  | grep Stepping | awk '{print $7}'`
  fi
  if [ -e t128$arga.out ]; then
    TIME3=`grep Time t128$arga.out | grep Stepping | awk '{print $7}'`
  fi
  if [ -e t128$argd.out ]; then
    TIME4=`grep Time t128$argd.out | grep Stepping | awk '{print $7}'`
  fi
  echo "$TIME1,$HOST1,$TIME2,$HOST2,$TIME3,$HOST3,$TIME4,$HOST4">>$OUTFILE
done
echo complete
