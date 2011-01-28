#!/bin/bash -f

# This script runs the a set of Validation Cases on a linux machine with
# a batch queuing system

export LD_LIBRARY_PATH=/shared/LIB64
export SVNROOT=`pwd`/../..
export FDS=$SVNROOT/FDS_Compilation/mpi_intel_linux_64/fds_mpi_intel_linux_64
export RUNFDSMPI=$SVNROOT/Utilities/Scripts/runfdsmpi.sh
export BASEDIR=`pwd`
export INDIR=FDS_Input_Files
# uncomment following line to stop all cases
#export STOPFDS=1

$RUNFDSMPI 5 $INDIR FM_SNL_04
$RUNFDSMPI 5 $INDIR FM_SNL_05
$RUNFDSMPI 5 $INDIR FM_SNL_21
