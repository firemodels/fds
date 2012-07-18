#!/bin/bash

# This script runs the FDS Verification Cases on a linux machine with
# a batch queuing system

export SVNROOT=`pwd`/..

# Set paths to FDS and FDS MPI executables
# If no argument is specfied, then run FDS release version.
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export FDSMPI=$SVNROOT/FDS_Compilation/mpi_intel_linux_64/fds_mpi_intel_linux_64

# Otherwise, if -d (debug) option is specified, then run FDS DB version.
while getopts 'dq:' OPTION
do
case $OPTION in
  d)
   export FDS=$SVNROOT/FDS_Compilation/intel_linux_64_db/fds_intel_linux_64_db
   export FDSMPI=$SVNROOT/FDS_Compilation/mpi_intel_linux_64_db/fds_mpi_intel_linux_64_db
   ;;
  q)
   queue="$OPTARG"
   ;;
esac
done

# VVVVVVVVVVVVVVVVV select which group of nodes to run on VVVVVVVV

# Set queue to submit cases to
# If no argument is specfied, then run cases on the default queue (batch).

   # blaze queue (default)
   export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
   export RUNFDSMPI=$SVNROOT/Utilities/Scripts/runfdsmpi.sh

# Otherwise, if -q (queue) option is specified, then run cases on the specified queue.
if [ "$queue" == "fire60s" ]
then
   # fire60s queue
   export RUNFDS=$SVNROOT/Utilities/Scripts/runfds6.sh
   export RUNFDSMPI=$SVNROOT/Utilities/Scripts/runfdsmpi6.sh
fi

if [ "$queue" == "fire70s" ]
then
   # fire70s queue
   export RUNFDS=$SVNROOT/Utilities/Scripts/runfds7.sh
   export RUNFDSMPI=$SVNROOT/Utilities/Scripts/runfdsmpi7.sh
fi

if [ "$queue" == "vis" ]
then
   # firevis/smokevis queue
   export RUNFDS=$SVNROOT/Utilities/Scripts/runfdsv.sh
   export RUNFDSMPI=$SVNROOT/Utilities/Scripts/runfdsmpiv.sh
fi

# ^^^^^^^^^^^^^^^^^ select nodes to run on ^^^^^^^^^^^^^^^^^^^^^^

export BASEDIR=`pwd`

# uncomment following line to stop all cases
# export STOPFDS=1

./FDS_Cases.sh
./FDS_MPI_Cases.sh

echo FDS cases submitted

