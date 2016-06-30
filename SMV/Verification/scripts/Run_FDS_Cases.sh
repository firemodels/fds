#!/bin/bash

# This script runs the FDS Verification Cases on a linux machine with
# a batch queuing system

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

if [ "$FDSNETWORK" == "infiniband" ] ; then
  IB=ib
fi

function usage {
echo "Run_FDS_Cases.sh [ -d -h -m max_iterations -o nthreads -q queue_name "
echo "                   -s -r resource_manager ]"
echo "Runs FDS verification suite"
echo ""
echo "Options"
echo "-d - use debug version of FDS"
echo "-h - display this message"
echo "-j - job prefix"
echo "-m max_iterations - stop FDS runs after a specifed number of iterations (delayed stop)"
echo "     example: an option of 10 would cause FDS to stop after 10 iterations"
echo "-o nthreads - run FDS with a specified number of threads [default: $nthreads]"
echo "-q queue_name - run cases using the queue queue_name"
echo "     default: batch"
echo "     other options: fire70s, vis"
echo "-r resource_manager - run cases using the resource manager"
echo "     default: PBS"
echo "     other options: SLURM"
echo "-s - stop FDS runs"
echo "-w time - walltime request for a batch job"
echo "     default: empty"
echo "     format for PBS: hh:mm:ss, format for SLURM: dd-hh:mm:ss"
exit
}

cd ../../..
export SVNROOT=`pwd`
cd $CURDIR

while getopts 'B:c:dD:hj:L:m:o:q:r:sw:' OPTION
do
case $OPTION in
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

IB=
if [ "$FDSNETWORK" == "infiniband" ]; then
  IB=ib
fi

export FDS=$SVNROOT/FDS/Build/${OPENMP}intel_$PLATFORM$DEBUG/fds_${OPENMP}intel_$PLATFORM$DEBUG
export FDSMPI=$SVNROOT/FDS/Build/mpi_intel_$PLATFORM$IB$DEBUG/fds_mpi_intel_$PLATFORM$IB$DEBUG
export QFDSSH="$SVNROOT/FDS/Utilities/Scripts/qfds.sh $RUNOPTION"

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

export QFDS="$QFDSSH $BACKGROUND $walltime -n $nthreads $JOBPREFIX -e $FDSMPI $QUEUE" 
cd ..
./FDS_Cases.sh
cd $CURDIR

echo FDS cases submitted
