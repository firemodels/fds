#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 96 -d $INDIR ethane_fireball.fds
$QFDS $DEBUG $QUEUE -p 96 -d $INDIR ethylene_fireball.fds
$QFDS $DEBUG $QUEUE -p 96 -d $INDIR isopentane_fireball.fds

echo FDS cases submitted
