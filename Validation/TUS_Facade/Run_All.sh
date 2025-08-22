#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 3 -d $INDIR TUS_facade_I-1-2.fds
$QFDS $DEBUG $QUEUE -p 3 -d $INDIR TUS_facade_1-1-2-3.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR TUS_facade_2-1-2-3.fds
$QFDS $DEBUG $QUEUE -p 27 -d $INDIR JIS_facade_5cm.fds
$QFDS $DEBUG $QUEUE -p 34 -d $INDIR JIS_facade_2cm.fds
$QFDS $DEBUG $QUEUE -p 90 -d $INDIR JIS_facade_1cm.fds
