#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR Bryant_034_kW.fds
$QFDS $DEBUG $QUEUE -d $INDIR Bryant_065_kW.fds
$QFDS $DEBUG $QUEUE -d $INDIR Bryant_096_kW.fds
$QFDS $DEBUG $QUEUE -d $INDIR Bryant_128_kW.fds
$QFDS $DEBUG $QUEUE -d $INDIR Bryant_160_kW.fds
$QFDS $DEBUG $QUEUE -d $INDIR Bryant_320_kW.fds
$QFDS $DEBUG $QUEUE -d $INDIR Bryant_511_kW.fds
 
echo FDS cases submitted
