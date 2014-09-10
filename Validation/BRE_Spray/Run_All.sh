#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR BRE_Spray_A_1.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Spray_A_2.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Spray_A_3.fds 
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Spray_A_4.fds 
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Spray_A_5.fds 
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Spray_A_6.fds 
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Spray_A_7.fds 
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Spray_A_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Spray_B_1.fds 
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Spray_B_2.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Spray_B_3.fds 
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Spray_B_4.fds 
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Spray_B_5.fds 
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Spray_B_6.fds 
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Spray_B_7.fds 
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Spray_B_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Spray_D_1.fds 
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Spray_D_2.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Spray_D_3.fds 
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Spray_D_4.fds 
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Spray_D_5.fds 
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Spray_D_6.fds 
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Spray_D_7.fds 
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Spray_D_8.fds
 
echo FDS cases submitted
