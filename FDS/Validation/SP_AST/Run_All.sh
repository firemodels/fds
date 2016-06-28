#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR SP_AST_Test_1.fds
$QFDS $DEBUG $QUEUE -d $INDIR SP_AST_Test_2.fds
$QFDS $DEBUG $QUEUE -d $INDIR SP_AST_Test_3.fds
$QFDS $DEBUG $QUEUE -d $INDIR SP_AST_Diesel_1p1.fds
$QFDS $DEBUG $QUEUE -d $INDIR SP_AST_Diesel_1p9.fds
$QFDS $DEBUG $QUEUE -d $INDIR SP_AST_Heptane_1p1.fds
$QFDS $DEBUG $QUEUE -d $INDIR SP_Room_A1.fds
$QFDS $DEBUG $QUEUE -d $INDIR SP_Room_A2.fds
$QFDS $DEBUG $QUEUE -d $INDIR SP_Room_A3.fds
$QFDS $DEBUG $QUEUE -d $INDIR SP_Room_A4.fds
$QFDS $DEBUG $QUEUE -d $INDIR SP_Room_A5.fds
$QFDS $DEBUG $QUEUE -d $INDIR SP_Room_B1.fds
$QFDS $DEBUG $QUEUE -d $INDIR SP_Room_B2.fds
$QFDS $DEBUG $QUEUE -d $INDIR SP_Room_C1.fds
$QFDS $DEBUG $QUEUE -d $INDIR SP_Room_C2.fds
$QFDS $DEBUG $QUEUE -d $INDIR SP_Room_C3.fds
$QFDS $DEBUG $QUEUE -d $INDIR SP_Room_D1.fds
$QFDS $DEBUG $QUEUE -d $INDIR SP_Room_D2.fds
$QFDS $DEBUG $QUEUE -d $INDIR plate_thermometer_25.fds
$QFDS $DEBUG $QUEUE -d $INDIR plate_thermometer_75.fds

echo FDS cases submitted
