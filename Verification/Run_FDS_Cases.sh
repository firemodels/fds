#!/bin/bash

# This script runs the FDS Verification Cases on a linux machine with
# a batch queuing system

queue=batch
cases=all
size=64
DEBUG=
IB=
nthreads=1
resource_manager=
export RUN_OPENMP

if [ "$FDSNETWORK" == "infiniband" ] ; then
IB=ib
fi

function usage {
echo "Run_FDS_Cases.sh [ -c cases -d -h -m max_iterations -o nthreads -q queue_name "
echo "                   -s -r resource_manager ]"
echo "Runs FDS verification suite"
echo ""
echo "Options"
echo "-c cases - select set of cases to run"
echo "     default: $cases"
echo "     other options: serial, mpi"
echo "-d - use debug version of FDS"
echo "-h - display this message"
echo "-m max_iterations - stop FDS runs after a specifed number of iterations (delayed stop)"
echo "     example: an option of 10 would cause FDS to stop after 10 iterations"
echo "-o nthreads - run OpenMP version of FDS with a specified number of threads [default: $nthreads]"
echo "-p size - platform size"
echo "     default: 64"
echo "     other options: 32"
echo "-q queue_name - run cases using the queue queue_name"
echo "     default: batch"
echo "     other options: fire60s, fire70s, vis"
echo "-r resource_manager - run cases using the resource manager"
echo "     default: PBS"
echo "     other options: SLURM"
echo "-s - stop FDS runs"
exit
}

export SVNROOT=`pwd`/..

while getopts 'c:dhm:o:p:q:r:s' OPTION
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
  m)
   export STOPFDSMAXITER="$OPTARG"
   ;;
  o)
   nthreads="$OPTARG"
   RUN_OPENMP=1
   ;;
  p)
   size="$OPTARG"
   ;;
  q)
   queue="$OPTARG"
   ;;
  r)
   resource_manager="$OPTARG"
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
if [ "$FDSNETWORK" == "infiniband" ]; then
  IB=ib
fi

export BACKGROUND=$SVNROOT/Utilities/background/intel_$PLATFORM2/background
export FDS=$SVNROOT/FDS_Compilation/${OPENMP}intel_$PLATFORM$DEBUG/fds_${OPENMP}intel_$PLATFORM$DEBUG
export FDSMPI=$SVNROOT/FDS_Compilation/mpi_intel_$PLATFORM$IB$DEBUG/fds_mpi_intel_$PLATFORM$IB$DEBUG

if [ "$resource_manager" == "SLURM" ]; then
   export RESOURCE_MANAGER="SLURM"
else
   export RESOURCE_MANAGER="PBS"
fi
if [ "$queue" != "" ]; then
   queue="-q $queue"
fi

export RUNFDS="$SVNROOT/Utilities/Scripts/runfds.sh $queue" 
export RUNFDSMPI="$SVNROOT/Utilities/Scripts/runfdsmpi.sh $queue"

if [ $RUN_OPENMP ]; then
  export RUNFDS="$SVNROOT/Utilities/Scripts/runfdsopenmp.sh $queue -n $nthreads"
  export RUNFDSBENCHMARK="$SVNROOT/Utilities/Scripts/runfdsopenmp.sh -b $queue"
else
  export RUNFDSBENCHMARK="$SVNROOT/Utilities/Scripts/runfds.sh $queue"
fi

export BASEDIR=`pwd`

# Run appropriate set of cases depending on user specified set (-c option)
case "$cases" in
  all)
    ./FDS_MPI_Cases.sh
    ./FDS_Cases.sh
    ;;
  serial)
    ./FDS_Cases.sh
    ;;
  mpi)
    ./FDS_MPI_Cases.sh
    ;;
esac

echo FDS cases submitted

