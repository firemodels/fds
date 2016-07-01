#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR Fleury_1t1_100_kW.fds
$QFDS $DEBUG $QUEUE -d $INDIR Fleury_1t1_150_kW.fds
$QFDS $DEBUG $QUEUE -d $INDIR Fleury_1t1_200_kW.fds
$QFDS $DEBUG $QUEUE -d $INDIR Fleury_1t1_250_kW.fds
$QFDS $DEBUG $QUEUE -d $INDIR Fleury_1t1_300_kW.fds
$QFDS $DEBUG $QUEUE -d $INDIR Fleury_2t1_100_kW.fds
$QFDS $DEBUG $QUEUE -d $INDIR Fleury_2t1_150_kW.fds
$QFDS $DEBUG $QUEUE -d $INDIR Fleury_2t1_200_kW.fds
$QFDS $DEBUG $QUEUE -d $INDIR Fleury_2t1_250_kW.fds
$QFDS $DEBUG $QUEUE -d $INDIR Fleury_2t1_300_kW.fds
$QFDS $DEBUG $QUEUE -d $INDIR Fleury_3t1_100_kW.fds
$QFDS $DEBUG $QUEUE -d $INDIR Fleury_3t1_150_kW.fds
$QFDS $DEBUG $QUEUE -d $INDIR Fleury_3t1_200_kW.fds
$QFDS $DEBUG $QUEUE -d $INDIR Fleury_3t1_250_kW.fds
$QFDS $DEBUG $QUEUE -d $INDIR Fleury_3t1_300_kW.fds

echo FDS cases submitted
