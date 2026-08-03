#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 55 -d $INDIR ethane_jet.fds
$QFDS $DEBUG $QUEUE -p 55 -d $INDIR ethylene_jet.fds
$QFDS $DEBUG $QUEUE -p 55 -d $INDIR isopentane_jet.fds

echo FDS cases submitted
