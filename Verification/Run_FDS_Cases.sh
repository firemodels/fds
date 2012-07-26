#!/bin/bash

# This script runs the FDS Verification Cases on a linux machine with
# a batch queuing system

queue=
cases=all

function usage {
echo "Run_FDS_Cases.sh [ -c cases -d -h -q queue_name -s ]"
echo "Runs FDS verification suite"
echo ""
echo "Options"
echo "-c cases - select set of cases to run"
echo "     default: $cases"
echo "     other options: serial, mpi"
echo "-d - use debug version of FDS"
echo "-h - display this message"
echo "-q queue_name - run cases using the queue queue_name"
echo "     default: batch"
echo "     other options: fire60s, fire70s, vis"
echo "-s - stop FDS runs"
exit
}

export SVNROOT=`pwd`/..

# Set paths to FDS and FDS MPI executables
# If no argument is specfied, then run FDS release version.
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export FDSMPI=$SVNROOT/FDS_Compilation/mpi_intel_linux_64/fds_mpi_intel_linux_64

# Otherwise, if -d (debug) option is specified, then run FDS DB version.
while getopts 'c:dhq:s' OPTION
do
case $OPTION in
  c)
   cases="$OPTARG"
   ;;
  d)
   export FDS=$SVNROOT/FDS_Compilation/intel_linux_64_db/fds_intel_linux_64_db
   export FDSMPI=$SVNROOT/FDS_Compilation/mpi_intel_linux_64_db/fds_mpi_intel_linux_64_db
   ;;
  h)
  usage;
  ;;
  q)
   queue="$OPTARG"
   ;;
  s)
   export STOPFDS=1
   ;;   
esac
done

if [ "$queue" != "" ]; then
   queue="-q $queue"
fi

export RUNFDS="$SVNROOT/Utilities/Scripts/runfds.sh $queue" 
export RUNFDSMPI="$SVNROOT/Utilities/Scripts/runfdsmpi.sh $queue"

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

