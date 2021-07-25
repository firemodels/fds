#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p  42 -d $INDIR jet1.fds
$QFDS $DEBUG $QUEUE -p  80 -d $INDIR jet2.fds
$QFDS $DEBUG $QUEUE -p 110 -d $INDIR jet3.fds

echo FDS cases submitted
