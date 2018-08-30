#/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 8  -d $INDIR UWO_test7_case1_180_32.fds
$QFDS $DEBUG $QUEUE -p 8  -d $INDIR UWO_test7_case1_180_64.fds
$QFDS $DEBUG $QUEUE -p 64 -d $INDIR UWO_test7_case1_180_128.fds
$QFDS $DEBUG $QUEUE -p 8  -d $INDIR UWO_test7_case1_270_32.fds
$QFDS $DEBUG $QUEUE -p 8  -d $INDIR UWO_test7_case1_270_64.fds
$QFDS $DEBUG $QUEUE -p 64 -d $INDIR UWO_test7_case1_270_128.fds

echo FDS cases submitted
