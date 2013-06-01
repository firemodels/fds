#!/bin/bash

# This script runs the FDS Verification Cases on a linux machine with
# a batch queuing system

queue=
cases=all
size=64
DEBUG=
OPENMP=
IB=
if [ "$FDSNETWORK" == "infiniband" ] ; then
IB=ib
fi

function usage {
echo "Run_FDS_Cases.sh [ -c cases -d -h -o -q queue_name -s ]"
echo "Runs FDS verification suite"
echo ""
echo "Options"
echo "-c cases - select set of cases to run"
echo "     default: $cases"
echo "     other options: serial, mpi"
echo "-d - use debug version of FDS"
echo "-h - display this message"
echo "-o - run OpenMP version of FDS"
echo "-p size - platform size"
echo "     default: 64"
echo "     other options: 32"
echo "-q queue_name - run cases using the queue queue_name"
echo "     default: batch"
echo "     other options: fire60s, fire70s, vis"
echo "-s - stop FDS runs"
exit
}

export SVNROOT=`pwd`/..

while getopts 'c:dhop:q:s' OPTION
do
case $OPTION in
  c)
   cases="$OPTARG"
   ;;
  d)
   DEBUG=_db
   ;;
  h)
  usage;
  ;;
  o)
  OPENMP=openmp_
  RUN_OPENMP=1
  ;;
  p)
  size="$OPTARG"
  ;;
  q)
   queue="$OPTARG"
   ;;
  s)
   export STOPFDS=1
   ;;   
esac
done

if [ "$size" != "32" ]; then
  size=64
fi
size=_$size

OS=`uname`
if [ "$OS" == "Darwin" ]; then
  PLATFORM=osx$size
  PLATFORM2=osx_32
else
  PLATFORM=linux$size
  PLATFORM2=linux_32
fi

IB=
if [ "$size" == "64" ]; then
if [ "$FDSNETWORK" == "infiniband" ]; then
  IB=ib
fi
fi

export BACKGROUND=$SVNROOT/Utilities/background/intel_$PLATFORM2/background
export FDS=$SVNROOT/FDS_Compilation/${OPENMP}intel_$PLATFORM$DEBUG/fds_${OPENMP}intel_$PLATFORM$DEBUG
export FDSMPI=$SVNROOT/FDS_Compilation/mpi_intel_$PLATFORM$IB$DEBUG/fds_mpi_intel_$PLATFORM$IB$DEBUG


if [ "$queue" != "" ]; then
   queue="-q $queue"
fi

export RUNFDS="$SVNROOT/Utilities/Scripts/runfds.sh $queue" 
export RUNFDSMPI="$SVNROOT/Utilities/Scripts/runfdsmpi.sh $queue"

if [ $RUN_OPENMP ]; then
  export RUNFDS="$SVNROOT/Utilities/Scripts/runfdsopenmp.sh $queue" 
fi

export BASEDIR=`pwd`

# Run appropriate set of cases depending on user specified set (-c option)
case "$cases" in
  all)
    ./FDS_Cases.sh
    ./FDS_MPI_Cases.sh
    ;;
  serial)
    ./FDS_Cases.sh
    ;;
  mpi)
    ./FDS_MPI_Cases.sh
    ;;
esac

echo FDS cases submitted

