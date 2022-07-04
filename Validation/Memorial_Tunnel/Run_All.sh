#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_606A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_612B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_615B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Cold_Flow_Series_1.fds

echo FDS cases submitted
