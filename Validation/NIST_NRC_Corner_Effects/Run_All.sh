#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 8 -d $INDIR corner_200_kW.fds
$QFDS $DEBUG $QUEUE -p 8 -d $INDIR corner_300_kW.fds
$QFDS $DEBUG $QUEUE -p 8 -d $INDIR corner_400_kW.fds

