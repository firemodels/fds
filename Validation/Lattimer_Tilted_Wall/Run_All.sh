#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR Lattimer_20_kW_0_degree_coarse.fds
$QFDS $DEBUG $QUEUE -d $INDIR Lattimer_20_kW_0_degree.fds
$QFDS $DEBUG -p 27 $QUEUE -d $INDIR Lattimer_20_kW_0_degree_fine.fds

$QFDS $DEBUG -p 4 $QUEUE -d $INDIR Lattimer_35_kW_00_degree.fds
$QFDS $DEBUG -p 4 $QUEUE -d $INDIR Lattimer_35_kW_67_degree.fds

$QFDS $DEBUG -p 4 $QUEUE -d $INDIR Lattimer_50_kW_00_degree.fds
$QFDS $DEBUG -p 4 $QUEUE -d $INDIR Lattimer_50_kW_11_degree.fds
$QFDS $DEBUG -p 4 $QUEUE -d $INDIR Lattimer_50_kW_23_degree.fds
$QFDS $DEBUG -p 4 $QUEUE -d $INDIR Lattimer_50_kW_35_degree.fds
$QFDS $DEBUG -p 4 $QUEUE -d $INDIR Lattimer_50_kW_49_degree.fds
$QFDS $DEBUG -p 4 $QUEUE -d $INDIR Lattimer_50_kW_67_degree.fds

$QFDS $DEBUG -p 4 $QUEUE -d $INDIR Lattimer_75_kW_00_degree.fds
$QFDS $DEBUG -p 4 $QUEUE -d $INDIR Lattimer_75_kW_67_degree.fds

$QFDS $DEBUG -p 16 $QUEUE -d $INDIR Lattimer_50_kW_00_degree_fine.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR Lattimer_50_kW_35_degree_fine.fds

echo FDS cases submitted
