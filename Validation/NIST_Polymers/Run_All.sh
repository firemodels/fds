#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR NIST_Gasification_Test_R3.fds
$QFDS $DEBUG $QUEUE -d $INDIR NIST_Gasification_Test_R4.fds
$QFDS $DEBUG $QUEUE -d $INDIR NIST_Gasification_Test_R5.fds
$QFDS $DEBUG $QUEUE -d $INDIR NIST_TGA_10K.fds
$QFDS $DEBUG $QUEUE -d $INDIR NIST_MCC_60K.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR NIST_Cone_Test_25kW.fds

echo FDS cases submitted
