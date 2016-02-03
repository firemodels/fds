#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02.fds
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_nodep.fds
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_5um.fds
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_10um.fds
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_ox.fds
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_ox_high.fds
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_ox_low.fds
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_5bin.fds
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_5bin_nodep.fds
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_5bin_ox.fds
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_5bin_ox_nodep.fds
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_5bin_ox_LR.fds
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_10bin.fds
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_10bin_ox.fds
$QFDS $DEBUG $QUEUE -p 4 -n 4 -d $INDIR  NIST_NRC_02_5bin_ox_HR.fds

echo FDS cases submitted
