#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_5bin_ox_2.fds
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_5bin_ox_5.fds
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_5bin_ox_10.fds
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_5bin_ox_20.fds

echo FDS cases submitted
