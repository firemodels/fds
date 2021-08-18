#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 5  $QUEUE -d $INDIR USN_Hawaii_Test_01.fds
$QFDS $DEBUG -p 5  $QUEUE -d $INDIR USN_Hawaii_Test_02.fds
$QFDS $DEBUG -p 5  $QUEUE -d $INDIR USN_Hawaii_Test_03.fds
$QFDS $DEBUG -p 5  $QUEUE -d $INDIR USN_Hawaii_Test_04.fds
$QFDS $DEBUG -p 5  $QUEUE -d $INDIR USN_Hawaii_Test_05.fds
$QFDS $DEBUG -p 5  $QUEUE -d $INDIR USN_Hawaii_Test_06.fds
$QFDS $DEBUG -p 5  $QUEUE -d $INDIR USN_Hawaii_Test_07.fds
$QFDS $DEBUG -p 5  $QUEUE -d $INDIR USN_Hawaii_Test_11.fds

$QFDS $DEBUG -p 5  $QUEUE -d $INDIR USN_Iceland_Test_01.fds
$QFDS $DEBUG -p 5  $QUEUE -d $INDIR USN_Iceland_Test_02.fds
$QFDS $DEBUG -p 5  $QUEUE -d $INDIR USN_Iceland_Test_03.fds
$QFDS $DEBUG -p 5  $QUEUE -d $INDIR USN_Iceland_Test_04.fds
$QFDS $DEBUG -p 5  $QUEUE -d $INDIR USN_Iceland_Test_05.fds
$QFDS $DEBUG -p 5  $QUEUE -d $INDIR USN_Iceland_Test_06.fds
$QFDS $DEBUG -p 5  $QUEUE -d $INDIR USN_Iceland_Test_07.fds
$QFDS $DEBUG -p 5  $QUEUE -d $INDIR USN_Iceland_Test_09.fds
$QFDS $DEBUG -p 5  $QUEUE -d $INDIR USN_Iceland_Test_10.fds
$QFDS $DEBUG -p 5  $QUEUE -d $INDIR USN_Iceland_Test_11.fds
$QFDS $DEBUG -p 5  $QUEUE -d $INDIR USN_Iceland_Test_12.fds
$QFDS $DEBUG -p 5  $QUEUE -d $INDIR USN_Iceland_Test_13.fds
$QFDS $DEBUG -p 5  $QUEUE -d $INDIR USN_Iceland_Test_14.fds
$QFDS $DEBUG -p 5  $QUEUE -d $INDIR USN_Iceland_Test_15.fds
$QFDS $DEBUG -p 5  $QUEUE -d $INDIR USN_Iceland_Test_17.fds
$QFDS $DEBUG -p 5  $QUEUE -d $INDIR USN_Iceland_Test_18.fds
$QFDS $DEBUG -p 5  $QUEUE -d $INDIR USN_Iceland_Test_19.fds
$QFDS $DEBUG -p 5  $QUEUE -d $INDIR USN_Iceland_Test_20.fds

