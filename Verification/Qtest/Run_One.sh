#!/bin/bash -f

# This script runs the a set of Validation Cases on a linux machine with
# a batch queuing system

export SVNROOT=`pwd`/../..
export FDSMPI=$SVNROOT/FDS_Compilation/mpi_intel_linux_64/fds_mpi_intel_linux_64
#export RUNFDSMPI=$SVNROOT/Utilities/Scripts/runfdsmpi.sh
export RUNFDSMPI=./runfdsmpi.sh
export BASEDIR=`pwd`
export INDIR=Current_Results
source ~/.bashrc_fds intel64

# uncomment following line to stop all cases
#export STOPFDS=1

/bin/sh -c "cp $BASEDIR/empty_box_16.fds $BASEDIR/$INDIR"

$RUNFDSMPI 16 $INDIR empty_box_16

echo FDS cases submitted
