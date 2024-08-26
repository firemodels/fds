#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

#$QFDS $DEBUG $QUEUE -p 864 -d $INDIR strong_scaling_test_864.fds
#$QFDS $DEBUG $QUEUE -p 864 -d $INDIR weak_scaling_test_864.fds

$QFDS $DEBUG $QUEUE -t -p 432 -d $INDIR weak_scaling_test_432.fds
$QFDS $DEBUG $QUEUE -t -p 288 -d $INDIR weak_scaling_test_288.fds
$QFDS $DEBUG $QUEUE -t -p 192 -d $INDIR weak_scaling_test_192.fds
$QFDS $DEBUG $QUEUE -t -p 128 -d $INDIR weak_scaling_test_128.fds
$QFDS $DEBUG $QUEUE -t -p 64  -d $INDIR weak_scaling_test_064.fds
$QFDS $DEBUG $QUEUE -t -p 32  -d $INDIR weak_scaling_test_032.fds
$QFDS $DEBUG $QUEUE -t -p 16  -d $INDIR weak_scaling_test_016.fds
$QFDS $DEBUG $QUEUE -t -p 8   -d $INDIR weak_scaling_test_008.fds
$QFDS $DEBUG $QUEUE -t -p 4   -d $INDIR weak_scaling_test_004.fds
$QFDS $DEBUG $QUEUE -t -p 2   -d $INDIR weak_scaling_test_002.fds
$QFDS $DEBUG $QUEUE -t -p 1   -d $INDIR weak_scaling_test_001.fds

$QFDS $DEBUG $QUEUE -t -p 432 -d $INDIR strong_scaling_test_432.fds
$QFDS $DEBUG $QUEUE -t -p 288 -d $INDIR strong_scaling_test_288.fds
$QFDS $DEBUG $QUEUE -t -p 192 -d $INDIR strong_scaling_test_192.fds
$QFDS $DEBUG $QUEUE -t -p 96  -d $INDIR strong_scaling_test_096.fds
$QFDS $DEBUG $QUEUE -t -p 64  -d $INDIR strong_scaling_test_064.fds
$QFDS $DEBUG $QUEUE -t -p 32  -d $INDIR strong_scaling_test_032.fds
$QFDS $DEBUG $QUEUE -t -p 8   -d $INDIR strong_scaling_test_008.fds
$QFDS $DEBUG $QUEUE -t -p 1   -d $INDIR strong_scaling_test_001.fds

