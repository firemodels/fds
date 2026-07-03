#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 41 -d $INDIR CERTEC_01_D3.fds
$QFDS $DEBUG $QUEUE -p 41 -d $INDIR CERTEC_04_D3.fds
$QFDS $DEBUG $QUEUE -p 41 -d $INDIR CERTEC_14_D4.fds
$QFDS $DEBUG $QUEUE -p 41 -d $INDIR CERTEC_10_D5.fds
$QFDS $DEBUG $QUEUE -p 41 -d $INDIR CERTEC_07_D6.fds

echo FDS cases submitted
