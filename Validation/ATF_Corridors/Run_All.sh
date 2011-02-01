#!/bin/bash

# This script runs the a set of Validation Cases on a linux machine with
# a batch queuing system

export SVNROOT=`pwd`/../..
export BASEDIR=`pwd`
export INDIR=Current_Results
source ../../FDS_Compilation/SET_MYFDSENV.sh intel64
#export FDSMPI=override FDSMPI defined in above script

# uncomment following line to stop all cases
#export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDSMPI 3 $INDIR ATF_Corridors_050_kW
$RUNFDSMPI 3 $INDIR ATF_Corridors_100_kW
$RUNFDSMPI 3 $INDIR ATF_Corridors_240_kW
$RUNFDSMPI 3 $INDIR ATF_Corridors_250_kW
$RUNFDSMPI 3 $INDIR ATF_Corridors_500_kW
$RUNFDSMPI 3 $INDIR ATF_Corridors_Mix_kW
