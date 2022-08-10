#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 1 $QUEUE -d $INDIR ribbed_channel_20.fds
$QFDS $DEBUG -p 1 $QUEUE -d $INDIR ribbed_channel_40.fds
$QFDS $DEBUG -p 24 $QUEUE -d $INDIR ribbed_channel_80.fds
$QFDS $DEBUG -p 192 $QUEUE -d $INDIR ribbed_channel_160.fds

$QFDS $DEBUG -p 1 $QUEUE -d $INDIR ribbed_channel_geom_20.fds
$QFDS $DEBUG -p 1 $QUEUE -d $INDIR ribbed_channel_geom_40.fds
$QFDS $DEBUG -p 24 $QUEUE -d $INDIR ribbed_channel_geom_80.fds
$QFDS $DEBUG -p 192 $QUEUE -d $INDIR ribbed_channel_geom_160.fds

echo FDS cases submitted
