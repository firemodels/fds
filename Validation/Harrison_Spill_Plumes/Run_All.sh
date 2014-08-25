#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS  -p 4 $QUEUE -d $INDIR SE4.fds
$QFDS  -p 4 $QUEUE -d $INDIR SE5.fds
$QFDS  -p 4 $QUEUE -d $INDIR SE6.fds
$QFDS  -p 4 $QUEUE -d $INDIR SE7.fds
$QFDS  -p 4 $QUEUE -d $INDIR SE8.fds
$QFDS  -p 4 $QUEUE -d $INDIR SE9.fds
$QFDS  -p 4 $QUEUE -d $INDIR SE10.fds
$QFDS  -p 4 $QUEUE -d $INDIR SE11.fds
$QFDS  -p 4 $QUEUE -d $INDIR SE12.fds
$QFDS  -p 4 $QUEUE -d $INDIR SE13.fds
$QFDS  -p 4 $QUEUE -d $INDIR SE14.fds
$QFDS  -p 4 $QUEUE -d $INDIR SE15.fds
$QFDS  -p 4 $QUEUE -d $INDIR SE16.fds
$QFDS  -p 4 $QUEUE -d $INDIR SE17.fds
$QFDS  -p 4 $QUEUE -d $INDIR SE18.fds
$QFDS  -p 4 $QUEUE -d $INDIR SE19.fds
$QFDS  -p 4 $QUEUE -d $INDIR SE20.fds
$QFDS  -p 4 $QUEUE -d $INDIR SE21.fds

echo FDS cases submitted
