#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 30 $QUEUE -d $INDIR Burro3.fds
$QFDS $DEBUG -p 30 $QUEUE -d $INDIR Burro7.fds
$QFDS $DEBUG -p 50 $QUEUE -d $INDIR Burro8.fds
$QFDS $DEBUG -p 30 $QUEUE -d $INDIR Burro9.fds
$QFDS $DEBUG -p 24 $QUEUE -d $INDIR Coyote3.fds
$QFDS $DEBUG -p 24 $QUEUE -d $INDIR Coyote5.fds
$QFDS $DEBUG -p 24 $QUEUE -d $INDIR Coyote6.fds
$QFDS $DEBUG -p 12 $QUEUE -d $INDIR Falcon1.fds
$QFDS $DEBUG -p 12 $QUEUE -d $INDIR Falcon3.fds
$QFDS $DEBUG -p 12 $QUEUE -d $INDIR Falcon4.fds
$QFDS $DEBUG -p 24 $QUEUE -d $INDIR MaplinSands27.fds
$QFDS $DEBUG -p 24 $QUEUE -d $INDIR MaplinSands34.fds
$QFDS $DEBUG -p 24 $QUEUE -d $INDIR MaplinSands35.fds

echo FDS cases submitted
