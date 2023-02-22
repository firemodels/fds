#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 27 $QUEUE -d $INDIR McCaffrey_14_kW_5.fds
$QFDS $DEBUG -p 27 $QUEUE -d $INDIR McCaffrey_22_kW_5.fds
$QFDS $DEBUG -p 27 $QUEUE -d $INDIR McCaffrey_33_kW_5.fds
$QFDS $DEBUG -p 27 $QUEUE -d $INDIR McCaffrey_45_kW_5.fds
$QFDS $DEBUG -p 27 $QUEUE -d $INDIR McCaffrey_57_kW_5.fds

$QFDS $DEBUG -p 27 $QUEUE -d $INDIR McCaffrey_14_kW_11.fds
$QFDS $DEBUG -p 27 $QUEUE -d $INDIR McCaffrey_22_kW_11.fds
$QFDS $DEBUG -p 27 $QUEUE -d $INDIR McCaffrey_33_kW_11.fds
$QFDS $DEBUG -p 27 $QUEUE -d $INDIR McCaffrey_45_kW_11.fds
$QFDS $DEBUG -p 27 $QUEUE -d $INDIR McCaffrey_57_kW_11.fds

$QFDS $DEBUG -p 125 $QUEUE -d $INDIR McCaffrey_14_kW_21.fds
$QFDS $DEBUG -p 125 $QUEUE -d $INDIR McCaffrey_22_kW_21.fds
$QFDS $DEBUG -p 125 $QUEUE -d $INDIR McCaffrey_33_kW_21.fds
$QFDS $DEBUG -p 125 $QUEUE -d $INDIR McCaffrey_45_kW_21.fds
$QFDS $DEBUG -p 125 $QUEUE -d $INDIR McCaffrey_57_kW_21.fds

$QFDS $DEBUG -p 158 $QUEUE -d $INDIR McCaffrey_14_kW_45.fds
$QFDS $DEBUG -p 158 $QUEUE -d $INDIR McCaffrey_22_kW_45.fds
$QFDS $DEBUG -p 158 $QUEUE -d $INDIR McCaffrey_33_kW_45.fds
$QFDS $DEBUG -p 158 $QUEUE -d $INDIR McCaffrey_45_kW_45.fds
$QFDS $DEBUG -p 158 $QUEUE -d $INDIR McCaffrey_57_kW_45.fds

echo FDS cases submitted
