#!/bin/bash

# This script runs the FDS Verification Cases on a linux machine with
# a batch queuing system

if [ ! -e .verification_script_dir ]; then
   echo "***error The Run_FDS_Cases.sh script must be run"
   echo "   in the directory Verification/script."
   echo "   script aborted."
   exit
fi

QUEUE=batch
DEBUG=
SINGLE=
nthreads=1
resource_manager=
walltime=
RUNOPTION=
CURDIR=`pwd`
QFDS_COUNT=/tmp/qfds_count_`whoami`
if [ "$BACKGROUND_PROG" == "" ]; then
  export BACKGROUND_PROG=background
fi
if [ "$BACKGROUND_DELAY" == "" ]; then
  BACKGROUND_DELAY=2
fi
if [ "$BACKGROUND_LOAD" == "" ]; then
  export BACKGROUND_LOAD=75
fi
REGULAR=1
BENCHMARK=1
OOPT=
POPT=
INTEL=
INTEL2=
GEOMCASES=
WAIT=

function usage {
echo "Run_FDS_Cases.sh [ -d -h -m max_iterations -o nthreads -q queue_name "
echo "                   -s -r resource_manager ]"
echo "Runs FDS verification suite"
echo ""
echo "Options"
echo "-b - run only benchmark cases"
echo "-d - use debug version of FDS"
echo "-g - run only geometry cases"
echo "-h - display this message"
echo "-j - job prefix"
echo "-J - use Intel MPI version of FDS"
echo "-m max_iterations - stop FDS runs after a specifed number of iterations (delayed stop)"
echo "     example: an option of 10 would cause FDS to stop after 10 iterations"
echo "-o nthreads - run FDS with a specified number of threads [default: $nthreads]"
echo "-q queue_name - run cases using the queue queue_name"
echo "     default: batch"
echo "     other options: fire70s, vis"
echo "-r resource_manager - run cases using the resource manager"
echo "     default: PBS"
echo "     other options: SLURM"
echo "-R - run only regular (non-benchmark) cases"
echo "-s - stop FDS runs"
echo "-w time - walltime request for a batch job"
echo "     default: empty"
echo "     format for PBS: hh:mm:ss, format for SLURM: dd-hh:mm:ss"
echo "-W - wait for cases to complete before returning"
exit
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

while getopts 'bB:c:dD:ghj:JL:m:o:q:r:RsS:w:W' OPTION
do
case $OPTION in
  b)
   BENCHMARK=1
   GEOMCASES=
   REGULAR=
   ;;
  d)
   DEBUG=_db
   SINGLE="1"
   ;;
  g)
   BENCHMARK=
   GEOMCASES=1
   REGULAR=
   ;;
  h)
   usage;
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
  o)
   nthreads="$OPTARG"
   ;;
  q)
   QUEUE="$OPTARG"
   ;;
  r)
   resource_manager="$OPTARG"
   ;;
  R)
   BENCHMARK=
   GEOMCASES=
   REGULAR=1
   ;;
  s)
   export STOPFDS=1
   ;;
  w)
   walltime="-w $OPTARG"
   ;;
  W)
   WAIT="1"
   ;;
esac
done

if [ "$JOBPREFIX" == "" ]; then
  JOBPREFIX=FB_
fi
export JOBPREFIX

size=_64

OS=`uname`
if [ "$OS" == "Darwin" ]; then
  PLATFORM=osx$size
else
  PLATFORM=linux$size
fi
if [ "$OOPT" != "" ]; then
  OOPT="-O $OOPT"
fi
if [ "$POPT" != "" ]; then
  POPT="-O $POPT"
fi

export FDS=$SVNROOT/fds/Build/${OPENMP}intel_$PLATFORM$DEBUG/fds_${OPENMP}intel_$PLATFORM$DEBUG
export FDSMPI=$SVNROOT/fds/Build/${INTEL}mpi_intel_$PLATFORM$DEBUG/fds_${INTEL}mpi_intel_$PLATFORM$DEBUG
export QFDSSH="$SVNROOT/fds/Utilities/Scripts/qfds.sh $RUNOPTION"

if [ "$resource_manager" == "SLURM" ]; then
   export RESOURCE_MANAGER="SLURM"
else
   export RESOURCE_MANAGER="PBS"
fi
if [ "$QUEUE" != "" ]; then
   if [ "$QUEUE" == "none" ]; then
      echo 0 > $QFDS_COUNT
   fi
   QUEUE="-q $QUEUE"
fi

export BASEDIR=`pwd`

export QFDS="$QFDSSH $walltime -n $nthreads $INTEL2 -e $FDSMPI $QUEUE $OOPT $POPT" 

cd ..
if [ "$BENCHMARK" == "1" ]; then
  if [ "$SINGLE" == "" ]; then
    ./FDS_Benchmark_Cases.sh
  else
    ./FDS_Benchmark_Cases_single.sh
  fi
  echo FDS benchmark cases submitted
fi

cd $CURDIR
cd ..
if [ "$REGULAR" == "1" ]; then
  ./FDS_Cases.sh
  echo FDS non-benchmark cases submitted
fi
cd $CURDIR
cd ..
if [ "$GEOMCASES" == "1" ]; then
  ./GEOM_Cases.sh
  echo FDS geometry cases submitted
fi
if [ "$WAIT" == "1" ]; then
  wait_cases_end
fi
cd $CURDIR
