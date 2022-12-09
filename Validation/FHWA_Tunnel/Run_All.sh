#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh
cp FDS_Input_Files/*txt Current_Results

$QFDS -p 12 $DEBUG $QUEUE -d $INDIR IFAB-07.fds
$QFDS -p 12 $DEBUG $QUEUE -d $INDIR IFAB-08.fds
$QFDS -p 12 $DEBUG $QUEUE -d $INDIR IFAB-09.fds
$QFDS -p 12 $DEBUG $QUEUE -d $INDIR IFAB-10.fds
$QFDS -p 12 $DEBUG $QUEUE -d $INDIR IFAB-11.fds
$QFDS -p 12 $DEBUG $QUEUE -d $INDIR IFAB-13.fds
$QFDS -p 12 $DEBUG $QUEUE -d $INDIR IFAB-14.fds
$QFDS -p 12 $DEBUG $QUEUE -d $INDIR IFAB-15.fds
$QFDS -p 12 $DEBUG $QUEUE -d $INDIR IFAB-19.fds
$QFDS -p 12 $DEBUG $QUEUE -d $INDIR IFAB-22.fds
$QFDS -p 12 $DEBUG $QUEUE -d $INDIR IFAB-24.fds

echo FDS cases submitted
