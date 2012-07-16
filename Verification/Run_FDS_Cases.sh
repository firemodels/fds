#!/bin/bash

# This script runs the FDS Verification Cases on a linux machine with
# a batch queuing system

export SVNROOT=`pwd`/..

# Set paths to FDS and FDS MPI executables
# If no argument is specfied, then run FDS release version.
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export FDSMPI=$SVNROOT/FDS_Compilation/mpi_intel_linux_64/fds_mpi_intel_linux_64

# Otherwise, if -d (debug) option is specified, then run FDS DB version.
while getopts 'd' OPTION
do
case $OPTION in
  d)
   export FDS=$SVNROOT/FDS_Compilation/intel_linux_64_db/fds_intel_linux_64_db
   export FDSMPI=$SVNROOT/FDS_Compilation/mpi_intel_linux_64_db/fds_mpi_intel_linux_64_db
   ;;
esac
done

# VVVVVVVVVVVVVVVVV select which group of nodes to run on VVVVVVVV

# new cluster
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export RUNFDSMPI=$SVNROOT/Utilities/Scripts/runfdsmpi.sh

# fire70s
#export RUNFDS=$SVNROOT/Utilities/Scripts/runfds7.sh
#export RUNFDSMPI=$SVNROOT/Utilities/Scripts/runfdsmpi7.sh

# firevis/smokevis
#export RUNFDS=$SVNROOT/Utilities/Scripts/runfdsv.sh
#export RUNFDSMPI=$SVNROOT/Utilities/Scripts/runfdsmpiv.sh

# ^^^^^^^^^^^^^^^^^ select nodes to run on ^^^^^^^^^^^^^^^^^^^^^^

export BASEDIR=`pwd`

# uncomment following line to stop all cases
# export STOPFDS=1

./FDS_Cases.sh
./FDS_MPI_Cases.sh

echo FDS cases submitted

