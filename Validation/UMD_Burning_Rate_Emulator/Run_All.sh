#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR BRE_Test_1_2D_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Test_2_2D_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Test_3_2D_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Test_4_2D_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Test_5_2D_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Test_6_2D_1mm.fds

$QFDS $DEBUG $QUEUE -d $INDIR BRE_Test_1_2D_p5mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Test_2_2D_p5mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Test_3_2D_p5mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Test_4_2D_p5mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Test_5_2D_p5mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Test_6_2D_p5mm.fds

echo FDS cases submitted
