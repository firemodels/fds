#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_501.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_502.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_605.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_606A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_607.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_608.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_610.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_611.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_612B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_615B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_617A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_618A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_621A.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_622B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_623B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_624B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Test_625B.fds
$QFDS -p 84 $DEBUG $QUEUE -d $INDIR Cold_Flow_Series_1.fds

echo FDS cases submitted
