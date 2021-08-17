#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR -p 24 acetone_1_m.fds
$QFDS $DEBUG $QUEUE -d $INDIR -p 24 ethanol_1_m.fds
$QFDS $DEBUG $QUEUE -d $INDIR -p 24 methanol_1_m.fds
$QFDS $DEBUG $QUEUE -d $INDIR -p 24 butane_1_m.fds
$QFDS $DEBUG $QUEUE -d $INDIR -p 24 benzene_1_m.fds
$QFDS $DEBUG $QUEUE -d $INDIR -p 24 heptane_1_m.fds
$QFDS $DEBUG $QUEUE -d $INDIR -p 1  ASTM_E2058_Water_Evap.fds
$QFDS $DEBUG $QUEUE -d $INDIR -p 96 VTT_heptane_1_m2.fds
$QFDS $DEBUG $QUEUE -d $INDIR -p 96 VTT_heptane_2_m2.fds

echo FDS cases submitted
