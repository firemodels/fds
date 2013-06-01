#!/bin/bash

# This script runs the a set of Validation Cases on a linux machine with
# a batch queuing system

export SVNROOT=`pwd`/../..
export QFDS=/usr/local/bin/qfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results
# qq="-q fire80s"
qq=

# uncomment following line to stop all cases
#export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$QFDS -r $qq -p 3 -d $INDIR ATF_Corridors_050_kW.fds
$QFDS -r $qq -p 3 -d $INDIR ATF_Corridors_100_kW.fds
$QFDS -r $qq -p 3 -d $INDIR ATF_Corridors_240_kW.fds
$QFDS -r $qq -p 3 -d $INDIR ATF_Corridors_250_kW.fds
$QFDS -r $qq -p 3 -d $INDIR ATF_Corridors_500_kW.fds
$QFDS -r $qq -p 3 -d $INDIR ATF_Corridors_Mix_kW.fds
