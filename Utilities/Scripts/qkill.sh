#!/bin/bash

# ---------------------------- usage ----------------------------------

function usage {
  echo "Usage: qkill.sh -b first job_id -e last job_id"
  echo ""
  echo "kill slurm jobs from first job_id to last job_id"
  echo ""
  echo " -b job_id - first job id"
  echo " -e job_id - last job id"
  echo " -h - show this message"
  echo ""
  echo "-b and -e parameters are required"
  exit
}

begin=
end=
abort=

while getopts 'b:e:h' OPTION
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
if [ "$abort" == "1" ]; then
  echo ""
  usage
  exit
fi


joblist=/tmp/joblist.$$
qstat -a | awk '{ print $1" "$2 }' | grep `whoami` | sort -u > $joblist
for job_id in `seq $begin $end`; do
   njobs=`grep $job_id $joblist | wc -l`
   if [ "$njobs" != "0" ]; then
     echo job $job_id canceled
     scancel $job_id
   fi
done
rm $joblist
