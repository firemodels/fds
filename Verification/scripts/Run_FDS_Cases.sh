#!/bin/bash

# This script runs the FDS Cases on a linux machine with
# a batch queuing system

if [ ! -e .verification_script_dir ]; then
   echo "***error The Run_FDS_Cases.sh script must be run"
   echo "   in the directory Verification/script."
   echo "   script aborted."
   exit
fi

if [ "$WINDIR" == "" ]; then
  QUEUE=batch
else
  QUEUE=terminal
fi
DEBUG=
SINGLE=
CURDIR=`pwd`
QFDS_COUNT=/tmp/qfds_count_`whoami`
if [ "$BACKGROUND_PROG" == "" ]; then
  export BACKGROUND_PROG=background
fi
if [ "$BACKGROUND_LOAD" == "" ]; then
  export BACKGROUND_LOAD=75
fi
# make Intel MPI the default on a linux system
INTEL=i
INTEL2="-I"
# make OpenMPI default on the Mac
if [ "`uname`" == "Darwin" ]; then
  INTEL=o
  INTEL2=
fi
WAIT=
CHECKCASES=
RESTART=
FDSEXEC=

function usage {
echo "Run_FDS_Cases.sh [ -d -h -m max_iterations -q queue_name -s "
echo "Runs FDS verification suite"
echo ""
echo "Options"
echo "-C - check that cases ran (used by firebot)"
echo "-d - use debug version of FDS"
echo "-e exe - override the full path of FDS used to run cases"
echo "-h - display this message"
echo "-j - job prefix"
echo "-J - use Intel MPI version of FDS"
echo "-m max_iterations - stop FDS runs after a specifed number of iterations (delayed stop)"
echo "     example: an option of 10 would cause FDS to stop after 10 iterations"
echo "-O - use OpenMPI version of FDS"
echo "-p - run picture cases"
echo "-q queue_name - run cases using the queue queue_name [default: batch]"
echo "-r - run restart test cases"
echo "-s - stop FDS runs"
echo "-W - wait for cases to complete before returning"
}

function get_full_path {
  filepath=$1

  if [[ $filepath == /* ]]; then
    full_filepath=$filepath
  else
    dir_filepath=$(dirname  "${filepath}")
    filename_filepath=$(basename  "${filepath}")
    curdir=`pwd`
    cd $dir_filepath
    full_filepath=`pwd`/$filename_filepath
    cd $curdir
  fi
}

wait_cases_end()
{
   if [[ "$QUEUE" == "none" ]]
   then
     while [[ `ps -u $USER -f | fgrep .fds | grep -v grep` != '' ]]; do
        JOBS_REMAINING=`ps -u $USER -f | fgrep .fds | grep -v grep | wc -l`
        echo "Waiting for ${JOBS_REMAINING} cases to complete."
        sleep 15
     done
   else
     while [[ `qstat -a | awk '{print $2 $4}' | grep $(whoami) | grep $JOBPREFIX` != '' ]]; do
        JOBS_REMAINING=`qstat -a | awk '{print $2 $4}' | grep $(whoami) | grep $JOBPREFIX | wc -l`
        echo "Waiting for ${JOBS_REMAINING} cases to complete."
        sleep 15
     done
   fi
}

cd ..
export SVNROOT=`pwd`/../..
cd $SVNROOT
export SVNROOT=`pwd`
cd $CURDIR
RUN_PICTURES=

while getopts 'Cdhj:e:Jm:Opq:rsW' OPTION
do
case $OPTION in
  C)
   CHECKCASES="1"
   ;;
  d)
   DEBUG="1"
   SINGLE="1"
   ;;
  e)
   FDSEXEC="-e $OPTARG"
   ;;
  h)
   usage;
   exit
   ;;
  j)
   JOBPREFIX="$OPTARG"
   ;;
  J)
   INTEL=i
   INTEL2="-I"
   ;;
  m)
   export STOPFDSMAXITER="$OPTARG"
   ;;
  O)
   INTEL=o
   INTEL2=
   ;;
  p)
   RUN_PICTURES=1
   RESTART=
   ;;
  q)
   QUEUE="$OPTARG"
   ;;
  r)
   RESTART=1
   RUN_PICTURES=
   ;;
  s)
   export STOPFDS=1
   ;;
  W)
   WAIT="1"
   ;;
  *)
   usage
   exit 1
   ;;
esac
done

if [ "$JOBPREFIX" == "" ]; then
  JOBPREFIX=FB_
fi
export JOBPREFIX

QFDSSH="$(pwd)/../../Utilities/Scripts/qfds.sh"

if [ "$QUEUE" != "" ]; then
   if [ "$QUEUE" == "none" ]; then
      echo 0 > $QFDS_COUNT
   fi
   QUEUE="-q $QUEUE"
fi

if [ "$DEBUG" != "" ]; then
   DEBUG="-T db"
fi

if [ "$CHECKCASES" == "1" ]; then
  QFDS="$(pwd)/Check_FDS_Cases.sh"
else
  QFDS="$QFDSSH $INTEL2 $QUEUE $DEBUG $FDSEXEC"
fi
export QFDS

cd $CURDIR
cd ..

if [[ "$RESTART" == "" ]] && [[ "$RUN_PICTURES" == "" ]]; then
   ./FDS_Cases.sh
   if [ "$CHECKCASES" == "" ]; then
      echo Cases in FDS_Cases.sh submitted
   fi
fi
if [[ "$RESTART" != "" ]] && [[ "$RUN_PICTURES" == "" ]]; then
   ./FDS_RESTART_Cases.sh 
   if [ "$CHECKCASES" == "" ]; then
      echo Cases in FDS_RESTART_Cases.sh submitted
   fi
fi
if [[ "$RUN_PICTURES" != "" ]]; then
  export RUNSMV=$QFDS
  ./FDS_Pictures.sh 
   if [ "$CHECKCASES" == "" ]; then
      echo Cases in FDS_Pictures.sh submitted
   fi
fi

cd $CURDIR
cd ..
if [ "$CHECKCASES" == "" ]; then
  if [ "$WAIT" == "1" ]; then
    wait_cases_end
  fi
fi
cd $CURDIR
