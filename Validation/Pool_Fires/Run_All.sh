#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh
#source $SVNROOT/Validation/Common_Run_All.sh -y -r SLURM -q serial -w 3-00:00:00

$QFDS $DEBUG $QUEUE -d $INDIR acetone_1_m.fds
$QFDS $DEBUG $QUEUE -d $INDIR ethanol_1_m.fds
$QFDS $DEBUG $QUEUE -d $INDIR methanol_1_m.fds
$QFDS $DEBUG $QUEUE -d $INDIR butane_1_m.fds
$QFDS $DEBUG $QUEUE -d $INDIR benzene_1_m.fds
$QFDS $DEBUG $QUEUE -d $INDIR heptane_1_m.fds
$QFDS $DEBUG $QUEUE -d $INDIR ASTM_E2058_Water_Evap.fds
$QFDS $DEBUG $QUEUE -d $INDIR -p 2 -n 2 -o 4 VTT_heptane_1_m2.fds
$QFDS $DEBUG $QUEUE -d $INDIR -p 2 -n 2 -o 4 VTT_heptane_2_m2.fds

echo FDS cases submitted
