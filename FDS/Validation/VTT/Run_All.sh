#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 5 $QUEUE -d $INDIR VTT_01.fds
$QFDS $DEBUG -p 5 $QUEUE -d $INDIR VTT_02.fds
$QFDS $DEBUG -p 5 $QUEUE -d $INDIR VTT_03.fds

echo FDS cases submitted

