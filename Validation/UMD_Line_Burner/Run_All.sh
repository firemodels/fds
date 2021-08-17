#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 24 $QUEUE -d $INDIR methane_dx_1p25cm.fds
$QFDS $DEBUG -p 24 $QUEUE -d $INDIR propane_dx_1p25cm.fds

$QFDS $DEBUG -p 164 $QUEUE -d $INDIR methane_dx_p625cm.fds
$QFDS $DEBUG -p 164 $QUEUE -d $INDIR propane_dx_p625cm.fds

$QFDS $DEBUG -p 204 $QUEUE -d $INDIR methane_dx_p3125cm.fds
$QFDS $DEBUG -p 204 $QUEUE -d $INDIR propane_dx_p3125cm.fds

echo FDS cases submitted
