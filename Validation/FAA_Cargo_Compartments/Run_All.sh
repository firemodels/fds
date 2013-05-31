#!/bin/bash -f

# This script runs the a set of Validation Cases on a linux machine with
# a batch queuing system

export SVNROOT=`pwd`/../..
export QFDS=/usr/local/bin/qfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results
# qq="-q fire80s"
qq=

# uncomment following line to stop all cases
export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$QFDS -r -p 5 $qq -d $INDIR  FAA_B707_Test_1.fds
$QFDS -r -p 5 $qq -d $INDIR  FAA_B707_Test_2.fds
$QFDS -r -p 5 $qq -d $INDIR  FAA_B707_Test_3.fds
$QFDS -r -p 5 $qq -d $INDIR  FAA_DC10.fds
