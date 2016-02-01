#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -d $INDIR $QUEUE 

$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_nodep
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_5um
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_10um
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_ox
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_ox_high
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_ox_low
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_5bin
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_5bin_nodep
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_5bin_ox
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_5bin_ox_nodep
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_5bin_ox_LR
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_10bin
$QFDS $DEBUG $QUEUE -d $INDIR  NIST_NRC_02_10bin_ox
$QFDS $DEBUG $QUEUE -p 4 -n 4 -d $INDIR  NIST_NRC_02_5bin_ox_HR

echo FDS cases submitted
