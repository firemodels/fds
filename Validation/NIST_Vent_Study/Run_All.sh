#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 8 -d $INDIR Test_1.fds
$QFDS $DEBUG $QUEUE -p 8 -d $INDIR Test_2.fds
$QFDS $DEBUG $QUEUE -p 8 -d $INDIR Test_3.fds
$QFDS $DEBUG $QUEUE -p 8 -d $INDIR Test_4.fds
$QFDS $DEBUG $QUEUE -p 8 -d $INDIR Test_5.fds
$QFDS $DEBUG $QUEUE -p 8 -d $INDIR Test_6.fds

