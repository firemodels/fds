#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS  $QUEUE -d $INDIR SP_AST_Test_1.fds
$QFDS  $QUEUE -d $INDIR SP_AST_Test_2.fds
$QFDS  $QUEUE -d $INDIR SP_AST_Test_3.fds
$QFDS  $QUEUE -d $INDIR SP_AST_Diesel_1p1.fds
$QFDS  $QUEUE -d $INDIR SP_AST_Diesel_1p9.fds
$QFDS  $QUEUE -d $INDIR SP_AST_Heptane_1p1.fds

echo FDS cases submitted
