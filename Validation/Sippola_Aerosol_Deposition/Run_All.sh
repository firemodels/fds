#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_01.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_02.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_03.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_04.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_05.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_06.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_07.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_08.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_09.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_10.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_11.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_12.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_13.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_14.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_15.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_16.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_17.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_18.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_19.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_20.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_21.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_22.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_23.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_24.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_25.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_26.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_27.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_28.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_29.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_30.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR Sippola_Test_31.fds

echo FDS cases submitted
