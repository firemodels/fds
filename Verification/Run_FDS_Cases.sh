#!/bin/bash

# This script runs the FDS Verification Cases on a linux machine with
# a batch queuing system

queue=batch

function usage {
echo "Run_FDS_Cases.sh [-d -h -q queue_name ]"
echo "Runs FDS verification suite"
echo ""
echo "Options"
echo "-d - use debug version of FDS"
echo "-h - display this message"
echo "-q queue_name - run cases using the queue queue_name"
echo "     default: $queue"
echo "     other options: fire60s, fire70s, vis"
exit
}

export SVNROOT=`pwd`/..

# Set paths to FDS and FDS MPI executables
# If no argument is specfied, then run FDS release version.
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export FDSMPI=$SVNROOT/FDS_Compilation/mpi_intel_linux_64/fds_mpi_intel_linux_64

# Otherwise, if -d (debug) option is specified, then run FDS DB version.
while getopts 'dhq:' OPTION
do
case $OPTION in
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
esac
done

export RUNFDS="$SVNROOT/Utilities/Scripts/runfds.sh -q $queue" 
export RUNFDSMPI="$SVNROOT/Utilities/Scripts/runfdsmpi.sh -q $queue"

export BASEDIR=`pwd`

# uncomment following line to stop all cases
# export STOPFDS=1

./FDS_Cases.sh
./FDS_MPI_Cases.sh

echo FDS cases submitted

