#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR NRCC_Facade_Win_1_05_MW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR NRCC_Facade_Win_1_06_MW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR NRCC_Facade_Win_1_08_MW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR NRCC_Facade_Win_2_05_MW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR NRCC_Facade_Win_2_06_MW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR NRCC_Facade_Win_2_08_MW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR NRCC_Facade_Win_2_10_MW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR NRCC_Facade_Win_3_05_MW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR NRCC_Facade_Win_3_06_MW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR NRCC_Facade_Win_3_08_MW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR NRCC_Facade_Win_3_10_MW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR NRCC_Facade_Win_4_05_MW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR NRCC_Facade_Win_4_06_MW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR NRCC_Facade_Win_4_08_MW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR NRCC_Facade_Win_4_10_MW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR NRCC_Facade_Win_5_05_MW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR NRCC_Facade_Win_5_06_MW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR NRCC_Facade_Win_5_08_MW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR NRCC_Facade_Win_5_10_MW.fds 

echo FDS cases submitted
