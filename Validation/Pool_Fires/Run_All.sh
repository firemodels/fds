#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR acetone_1_m.fds
$QFDS $DEBUG $QUEUE -d $INDIR ethanol_1_m.fds
$QFDS $DEBUG $QUEUE -d $INDIR methanol_1_m.fds
$QFDS $DEBUG $QUEUE -d $INDIR butane_1_m.fds
$QFDS $DEBUG $QUEUE -d $INDIR benzene_1_m.fds
$QFDS $DEBUG $QUEUE -d $INDIR heptane_1_m.fds

echo FDS cases submitted
