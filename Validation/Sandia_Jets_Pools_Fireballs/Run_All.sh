#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 144 -d $INDIR ethane_pool.fds
$QFDS $DEBUG $QUEUE -p 144 -d $INDIR ethylene_pool.fds
$QFDS $DEBUG $QUEUE -p 144 -d $INDIR isopentane_pool.fds
$QFDS $DEBUG $QUEUE -p 144 -d $INDIR propane_pool.fds

echo FDS cases submitted
