#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 2 $QUEUE -d $INDIR  NBS_100A.fds
$QFDS $DEBUG -p 2 $QUEUE -d $INDIR  NBS_100O.fds
$QFDS $DEBUG -p 3 $QUEUE -d $INDIR  NBS_100Z.fds
