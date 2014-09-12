#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR McCaffrey_14_kW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR McCaffrey_22_kW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR McCaffrey_33_kW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR McCaffrey_45_kW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR McCaffrey_57_kW.fds 

$QFDS $DEBUG $QUEUE -d $INDIR McCaffrey_14_kW_coarse.fds
$QFDS $DEBUG $QUEUE -d $INDIR McCaffrey_22_kW_coarse.fds
$QFDS $DEBUG $QUEUE -d $INDIR McCaffrey_33_kW_coarse.fds
$QFDS $DEBUG $QUEUE -d $INDIR McCaffrey_45_kW_coarse.fds
$QFDS $DEBUG $QUEUE -d $INDIR McCaffrey_57_kW_coarse.fds

$QFDS $DEBUG -p 4 $QUEUE -d $INDIR McCaffrey_14_kW_fine.fds
$QFDS $DEBUG -p 4 $QUEUE -d $INDIR McCaffrey_22_kW_fine.fds
$QFDS $DEBUG -p 4 $QUEUE -d $INDIR McCaffrey_33_kW_fine.fds
$QFDS $DEBUG -p 4 $QUEUE -d $INDIR McCaffrey_45_kW_fine.fds
$QFDS $DEBUG -p 4 $QUEUE -d $INDIR McCaffrey_57_kW_fine.fds

#$QFDS $DEBUG $QUEUE -d $INDIR McCaffrey_14_kW_5.fds 
#$QFDS $DEBUG $QUEUE -d $INDIR McCaffrey_22_kW_5.fds 
#$QFDS $DEBUG $QUEUE -d $INDIR McCaffrey_33_kW_5.fds 
#$QFDS $DEBUG $QUEUE -d $INDIR McCaffrey_45_kW_5.fds 
#$QFDS $DEBUG $QUEUE -d $INDIR McCaffrey_57_kW_5.fds

#$QFDS $DEBUG -p 12 $QUEUE -d $INDIR McCaffrey_14_kW_10.fds
#$QFDS $DEBUG -p 12 $QUEUE -d $INDIR McCaffrey_22_kW_10.fds
#$QFDS $DEBUG -p 12 $QUEUE -d $INDIR McCaffrey_33_kW_10.fds
#$QFDS $DEBUG -p 12 $QUEUE -d $INDIR McCaffrey_45_kW_10.fds
#$QFDS $DEBUG -p 12 $QUEUE -d $INDIR McCaffrey_57_kW_10.fds

echo FDS cases submitted
