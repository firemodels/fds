#!/bin/bash -f

# This script runs the a set of Validation Cases on a linux machine with
# a batch queuing system

export SVNROOT=`pwd`/../..
export FDS=$SVNROOT/FDS_Compilation/mpi_intel_linux_64/fds_mpi_intel_linux_64
export RUNFDSMPI=$SVNROOT/Utilities/Scripts/runfdsmpi.sh
export BASEDIR=`pwd`
export INDIR=Current_Results

# uncomment following line to stop all cases
# export STOPFDS=1

cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR

$RUNFDSMPI 3 $INDIR ATF_Corridors_050_kW
$RUNFDSMPI 3 $INDIR ATF_Corridors_100_kW
$RUNFDSMPI 3 $INDIR ATF_Corridors_240_kW
$RUNFDSMPI 3 $INDIR ATF_Corridors_250_kW
$RUNFDSMPI 3 $INDIR ATF_Corridors_500_kW
$RUNFDSMPI 3 $INDIR ATF_Corridors_Mix_kW


