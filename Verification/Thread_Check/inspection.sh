#!/bin/bash
CURDIR=`pwd`
cd ../../..
GITROOT=`pwd`
cd $CURDIR
RESULT_DIR=$CURDIR/inspect_results
PROCESSES=1
QFDS="$GITROOT/fds/Utilities/Scripts/qfds.sh"
QUEUE=

function usage {
  echo "Usage: inspect_openmp.sh [-v] casename.fds"
  echo ""
  echo " -h display this message"
  echo " -v   - list command that will be used to thread check"
  echo " -d   - select directory to store Inspector results"
  echo " -p   - select number of MPI processes"
  echo "input_file - input file"
  echo ""
  exit
}

if [ $# -lt 1 ]
then
  usage
fi

showinput=
while getopts 'hd:vp:q:' OPTION
do
case $OPTION  in
  h)
   usage;
   ;;
  d)
   RESULT_DIR=$OPTARG
   ;;
  v)
   showinput=1
   ;;
  p)
   PROCESSES=$OPTARG
   ;;
  q)
   QUEUE=$OPTARG
esac
done
shift $(($OPTIND-1))
case=$1

# Perform OpenMP thread checking (detect deadlocks and data races)

cd $CURDIR
$QFDS -p $PROCESSES -o 4 -x $RESULT_DIR $case
sleep 5


while [[ `qstat -a | awk '{print $2" "$4}' | grep $(whoami) | grep '${case%.fds}'` != '' ]]; do
        echo "Waiting for case to complete." >> $CURDIR/monitor.txt
        sleep 240
done

if [[ `grep -i -E 'warning|remark|problem|error' $CURDIR/${case%.fds}.err | grep -v '0 new problem(s) found' | grep -v 'Warning: One or more threads in the application accessed the stack of another thread'` == "" ]]
   then
      echo "Inspector Clean" >> $CURDIR/monitor.txt
   else
      echo "Inspector found errors, use Inspector to analyze results." >> $CURDIR/monitor.txt
    
   fi
