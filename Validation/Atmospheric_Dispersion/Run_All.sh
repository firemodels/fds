#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 12 -d $INDIR plume_rise_1.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR plume_rise_2.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR plume_rise_3.fds

echo FDS cases submitted
