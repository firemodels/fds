#/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 8   -n 8  -d $INDIR UWO_Test7_Case1_180_5p00mm.fds
$QFDS $DEBUG $QUEUE -p 64  -n 8  -d $INDIR UWO_Test7_Case1_180_2p50mm.fds
$QFDS $DEBUG $QUEUE -p 120 -n 8  -d $INDIR UWO_Test7_Case1_180_1p25mm.fds
$QFDS $DEBUG $QUEUE -p 8   -n 8  -d $INDIR UWO_Test7_Case1_270_5p00mm.fds
$QFDS $DEBUG $QUEUE -p 64  -n 8  -d $INDIR UWO_Test7_Case1_270_2p50mm.fds
$QFDS $DEBUG $QUEUE -p 120 -n 8  -d $INDIR UWO_Test7_Case1_270_1p25mm.fds

echo FDS cases submitted
