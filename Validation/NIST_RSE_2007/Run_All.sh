#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 2 -d $INDIR RSE_01.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR RSE_02.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR RSE_03.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR RSE_04.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR RSE_05.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR RSE_06.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR RSE_07.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR RSE_10.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR RSE_11.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR RSE_12.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR RSE_15.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR RSE_16.fds

echo FDS cases submitted
