#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 126 -d $INDIR TAMU_T16.fds
$QFDS $DEBUG $QUEUE -p 126 -d $INDIR TAMU_T17.fds
$QFDS $DEBUG $QUEUE -p 126 -d $INDIR TAMU_T18.fds
$QFDS $DEBUG $QUEUE -p 126 -d $INDIR TAMU_T19a.fds
$QFDS $DEBUG $QUEUE -p 126 -d $INDIR TAMU_T19b.fds
$QFDS $DEBUG $QUEUE -p 126 -d $INDIR TAMU_T19c.fds
$QFDS $DEBUG $QUEUE -p 126 -d $INDIR TAMU_T19d.fds
$QFDS $DEBUG $QUEUE -p 126 -d $INDIR TAMU_T19e.fds
$QFDS $DEBUG $QUEUE -p 126 -d $INDIR TAMU_T19f.fds

echo FDS cases submitted
