#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR Lattimer_20_kW_0_degree_coarse.fds
$QFDS $DEBUG $QUEUE -d $INDIR Lattimer_20_kW_0_degree.fds
$QFDS $DEBUG -p 27 $QUEUE -d $INDIR Lattimer_20_kW_0_degree_fine.fds

$QFDS $DEBUG -p 16 $QUEUE -d $INDIR Lattimer_50_kW_00_degree_fine.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR Lattimer_50_kW_35_degree_fine.fds

echo FDS cases submitted
