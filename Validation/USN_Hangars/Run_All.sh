#!/bin/bash

# This script runs the a set of Validation Cases on a linux machine with
# a batch queuing system

export SVNROOT=`pwd`/../..
export QFDS=/usr/local/bin/qfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results
# qq="-q fire80s"
qq=

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$QFDS -p 5 -r $qq -d $INDIR USN_Hawaii_Test_01.fds
$QFDS -p 5 -r $qq -d $INDIR USN_Hawaii_Test_02.fds
$QFDS -p 5 -r $qq -d $INDIR USN_Hawaii_Test_03.fds
$QFDS -p 5 -r $qq -d $INDIR USN_Hawaii_Test_04.fds
$QFDS -p 5 -r $qq -d $INDIR USN_Hawaii_Test_05.fds
$QFDS -p 5 -r $qq -d $INDIR USN_Hawaii_Test_06.fds
$QFDS -p 5 -r $qq -d $INDIR USN_Hawaii_Test_07.fds
$QFDS -p 5 -r $qq -d $INDIR USN_Hawaii_Test_11.fds

$QFDS -p 5 -r $qq -d $INDIR USN_Iceland_Test_01.fds
$QFDS -p 5 -r $qq -d $INDIR USN_Iceland_Test_02.fds
$QFDS -p 5 -r $qq -d $INDIR USN_Iceland_Test_03.fds
$QFDS -p 5 -r $qq -d $INDIR USN_Iceland_Test_04.fds
$QFDS -p 5 -r $qq -d $INDIR USN_Iceland_Test_05.fds
$QFDS -p 5 -r $qq -d $INDIR USN_Iceland_Test_06.fds
$QFDS -p 5 -r $qq -d $INDIR USN_Iceland_Test_07.fds
$QFDS -p 5 -r $qq -d $INDIR USN_Iceland_Test_09.fds
$QFDS -p 5 -r $qq -d $INDIR USN_Iceland_Test_10.fds
$QFDS -p 5 -r $qq -d $INDIR USN_Iceland_Test_11.fds
$QFDS -p 5 -r $qq -d $INDIR USN_Iceland_Test_12.fds
$QFDS -p 5 -r $qq -d $INDIR USN_Iceland_Test_13.fds
$QFDS -p 5 -r $qq -d $INDIR USN_Iceland_Test_14.fds
$QFDS -p 5 -r $qq -d $INDIR USN_Iceland_Test_15.fds
$QFDS -p 5 -r $qq -d $INDIR USN_Iceland_Test_17.fds
$QFDS -p 5 -r $qq -d $INDIR USN_Iceland_Test_18.fds
$QFDS -p 5 -r $qq -d $INDIR USN_Iceland_Test_19.fds
$QFDS -p 5 -r $qq -d $INDIR USN_Iceland_Test_20.fds

