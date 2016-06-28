#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR Sippola_Test_01.fds
$QFDS $DEBUG $QUEUE -d $INDIR Sippola_Test_02.fds
$QFDS $DEBUG $QUEUE -d $INDIR Sippola_Test_03.fds
$QFDS $DEBUG $QUEUE -d $INDIR Sippola_Test_04.fds
$QFDS $DEBUG $QUEUE -d $INDIR Sippola_Test_05.fds
$QFDS $DEBUG $QUEUE -d $INDIR Sippola_Test_06.fds
$QFDS $DEBUG $QUEUE -d $INDIR Sippola_Test_07.fds
$QFDS $DEBUG $QUEUE -d $INDIR Sippola_Test_08.fds
$QFDS $DEBUG $QUEUE -d $INDIR Sippola_Test_09.fds
$QFDS $DEBUG $QUEUE -d $INDIR Sippola_Test_10.fds
$QFDS $DEBUG $QUEUE -d $INDIR Sippola_Test_11.fds
$QFDS $DEBUG $QUEUE -d $INDIR Sippola_Test_12.fds
$QFDS $DEBUG $QUEUE -d $INDIR Sippola_Test_13.fds
$QFDS $DEBUG $QUEUE -d $INDIR Sippola_Test_14.fds
$QFDS $DEBUG $QUEUE -d $INDIR Sippola_Test_15.fds
$QFDS $DEBUG $QUEUE -d $INDIR Sippola_Test_16.fds

echo FDS cases submitted
