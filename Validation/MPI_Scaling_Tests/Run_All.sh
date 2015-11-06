#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 432 -n 12 -d $INDIR weak_scaling_test_432.fds
$QFDS $DEBUG $QUEUE -p 288 -n 12 -d $INDIR weak_scaling_test_288.fds
$QFDS $DEBUG $QUEUE -p 192 -n 12 -d $INDIR weak_scaling_test_192.fds
$QFDS $DEBUG $QUEUE -p 128 -n 12 -d $INDIR weak_scaling_test_128.fds
$QFDS $DEBUG $QUEUE -p 64  -n 12 -d $INDIR weak_scaling_test_064.fds
$QFDS $DEBUG $QUEUE -p 32  -n 12 -d $INDIR weak_scaling_test_032.fds
$QFDS $DEBUG $QUEUE -p 16  -n 12 -d $INDIR weak_scaling_test_016.fds
$QFDS $DEBUG $QUEUE -p 8   -n 12 -d $INDIR weak_scaling_test_008.fds
$QFDS $DEBUG $QUEUE -p 4   -n 12 -d $INDIR weak_scaling_test_004.fds
$QFDS $DEBUG $QUEUE -p 2   -n 12 -d $INDIR weak_scaling_test_002.fds
$QFDS $DEBUG $QUEUE -p 1   -n 12 -d $INDIR weak_scaling_test_001.fds

$QFDS $DEBUG $QUEUE -p 432 -n 12 -d $INDIR strong_scaling_test_432.fds
$QFDS $DEBUG $QUEUE -p 288 -n 12 -d $INDIR strong_scaling_test_288.fds
$QFDS $DEBUG $QUEUE -p 192 -n 12 -d $INDIR strong_scaling_test_192.fds
$QFDS $DEBUG $QUEUE -p 96  -n 12 -d $INDIR strong_scaling_test_096.fds
$QFDS $DEBUG $QUEUE -p 64  -n 12 -d $INDIR strong_scaling_test_064.fds
$QFDS $DEBUG $QUEUE -p 32  -n 12 -d $INDIR strong_scaling_test_032.fds
$QFDS $DEBUG $QUEUE -p 8   -n 12 -d $INDIR strong_scaling_test_008.fds
$QFDS $DEBUG $QUEUE -p 1   -n 12 -d $INDIR strong_scaling_test_001.fds

