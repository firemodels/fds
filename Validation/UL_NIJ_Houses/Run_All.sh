#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Single_Story_Gas_1.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Single_Story_Gas_2.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Single_Story_Gas_5.fds
$QFDS $DEBUG $QUEUE -p 40 -d $INDIR Two_Story_Gas_1.fds
$QFDS $DEBUG $QUEUE -p 40 -d $INDIR Two_Story_Gas_4.fds
$QFDS $DEBUG $QUEUE -p 40 -d $INDIR Two_Story_Gas_6.fds

echo FDS cases submitted
