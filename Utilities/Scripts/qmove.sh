#!/bin/bash

# ---------------------------- usage ----------------------------------

function usage {
  echo "Usage: qmove.sh -b first job_id -e last job_id [-q queue]"
  echo ""
  echo "move slurm jobs from first job_id to last job_id to the queue specified by -q"
  echo ""
  echo " -b job_id - first job id"
  echo " -e job_id - last job id"
  echo " -h - show this message"
  echo " -q queue - move jobs to queue"
  echo ""
  echo "-b, -e and -q parameters are required"
  exit
}

queue=
begin=
end=
abort=

while getopts 'b:e:hq:' OPTION
do
case $OPTION  in
  b) 
   begin="$OPTARG"
   ;;
  e) 
   end="$OPTARG"
   ;;
  h) 
   usage
   ;;
  q) 
   queue="$OPTARG"
   ;;
esac
done
shift $(($OPTIND-1))

if [ "$begin" == "" ]; then
  echo "***error: -b parameter not specified"
  abort=1
fi
if [ "$end" == "" ]; then
  echo "***error: -e parameter not specified"
  abort=1
fi
if [ "$queue" == "" ]; then
  echo "***error: -q parameter not specified"
  abort=1
fi
if [ "$abort" == "1" ]; then
  echo ""
  usage
  exit
fi


joblist=/tmp/joblist.$$
qstat -a | awk '{ print $1" "$2" "$10}' | grep `whoami` | grep 'Q$' | sort -u > $joblist
for job_id in `seq $begin $end`; do
   njobs=`grep $job_id $joblist | wc -l`
   if [ "$njobs" != "0" ]; then
     echo job $job_id moved to queue $queue
     scontrol update jobid=$job_id partition=$queue
   fi
done
rm $joblist
