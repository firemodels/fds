#!/bin/bash -f

# This script runs the a set of Validation Cases on a linux machine with
# a batch queuing system

export SVNROOT=`pwd`/../..
export FDSMPI=$SVNROOT/FDS_Compilation/mpi_intel_linux_64/fds_mpi_intel_linux_64
export RUNFDSMPI=$SVNROOT/Utilities/Scripts/runfdsmpi.sh
export BASEDIR=`pwd`
export INDIR=../Visualization2
source ~/.bashrc_fds intel64

# uncomment following line to stop all cases
#export STOPFDS=1

#$RUNFDSMPI 16 $INDIR vis_mtest1
#$RUNFDSMPI 16 $INDIR vis_mtest2
#$RUNFDSMPI 16 $INDIR vis_mtest3
$RUNFDSMPI 128 $INDIR mplume5c128
$RUNFDSMPI 8 $INDIR mplume5c8
