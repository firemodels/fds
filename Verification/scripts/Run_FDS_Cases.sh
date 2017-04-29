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
IB=
nthreads=1
resource_manager=
walltime=
RUNOPTION=
CURDIR=`pwd`
BACKGROUND=
BACKGROUND_DELAY=2
BACKGROUND_LOAD=75
JOBPREFIX=
REGULAR=1
BENCHMARK=1
OOPT=
POPT=

if [ "$FDSNETWORK" == "infiniband" ] ; then
  IB=ib
fi

function usage {
echo "Run_FDS_Cases.sh [ -d -h -m max_iterations -o nthreads -q queue_name "
echo "                   -s -r resource_manager ]"
echo "Runs FDS verification suite"
echo ""
echo "Options"
echo "-b - run only benchmark cases"
echo "-d - use debug version of FDS"
echo "-h - display this message"
echo "-j - job prefix"
echo "-m max_iterations - stop FDS runs after a specifed number of iterations (delayed stop)"
echo "     example: an option of 10 would cause FDS to stop after 10 iterations"
echo "-o nthreads - run FDS with a specified number of threads [default: $nthreads]"
echo "-O - pass through -O option to qfds.sh"
echo "-P - pass through -P option to qfds.sh"
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
exit
}

cd ..
export SVNROOT=`pwd`/../..
cd $SVNROOT
export SVNROOT=`pwd`
cd $CURDIR

while getopts 'bB:c:dD:hj:L:m:o:O:P:q:r:Rsw:' OPTION
do
case $OPTION in
  b)
   BENCHMARK=1
   REGULAR=
   ;;
  R)
   BENCHMARK=
   REGULAR=1
   ;;
  d)
   DEBUG=_db
   ;;
  B)
   BACKGROUND="$OPTARG"
   ;;
  D)
   BACKGROUND_DELAY="$OPTARG"
   ;;
  h)
   usage;
   ;;
  j)
   JOBPREFIX="-j $OPTARG"
   ;;
  L)
   BACKGROUND_LOAD="$OPTARG"
   ;;
  m)
   export STOPFDSMAXITER="$OPTARG"
   ;;
  o)
   nthreads="$OPTARG"
   ;;
  O)
   OOPT="$OPTARG"
   ;;
  P)
   POPT="$OPTARG"
   ;;
  q)
   QUEUE="$OPTARG"
   ;;
  r)
   resource_manager="$OPTARG"
   ;;
  s)
   export STOPFDS=1
   ;;
  w)
   walltime="-w $OPTARG"
   ;;
esac
done

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

IB=
if [ "$FDSNETWORK" == "infiniband" ]; then
  IB=ib
fi

export FDS=$SVNROOT/fds/Build/${OPENMP}intel_$PLATFORM$DEBUG/fds_${OPENMP}intel_$PLATFORM$DEBUG
export FDSMPI=$SVNROOT/fds/Build/mpi_intel_$PLATFORM$IB$DEBUG/fds_mpi_intel_$PLATFORM$IB$DEBUG
export QFDSSH="$SVNROOT/fds/Utilities/Scripts/qfds.sh $RUNOPTION"

if [ "$resource_manager" == "SLURM" ]; then
   export RESOURCE_MANAGER="SLURM"
else
   export RESOURCE_MANAGER="PBS"
fi
if [ "$QUEUE" != "" ]; then
   if [ "$QUEUE" == "none" ]; then
      if [ "$BACKGROUND" == "" ]; then
         BACKGROUND=background
      fi
      BACKGROUND="-B $BACKGROUND"
      export BACKGROUND_DELAY
      export BACKGROUND_LOAD
      JOBPREFIX=
   fi
   QUEUE="-q $QUEUE"
fi

export BASEDIR=`pwd`

export QFDS="$QFDSSH $BACKGROUND $walltime -n $nthreads $JOBPREFIX -e $FDSMPI $QUEUE $OOPT $POPT" 

cd ..
if [ "$BENCHMARK" == "1" ]; then
  ./FDS_Benchmark_Cases.sh
  echo FDS benchmark cases submitted
fi

cd $CURDIR
cd ..
if [ "$REGULAR" == "1" ]; then
  ./FDS_Cases.sh
  echo FDS non-benchmark cases submitted
fi
cd $CURDIR
