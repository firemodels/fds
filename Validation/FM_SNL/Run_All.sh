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
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$QFDS -r -p 5 $qq -d $INDIR FM_SNL_01.fds
$QFDS -r -p 5 $qq -d $INDIR FM_SNL_02.fds
$QFDS -r -p 5 $qq -d $INDIR FM_SNL_03.fds
$QFDS -r -p 5 $qq -d $INDIR FM_SNL_04.fds
$QFDS -r -p 5 $qq -d $INDIR FM_SNL_05.fds
$QFDS -r -p 5 $qq -d $INDIR FM_SNL_06.fds
$QFDS -r -p 5 $qq -d $INDIR FM_SNL_07.fds
$QFDS -r -p 5 $qq -d $INDIR FM_SNL_08.fds
$QFDS -r -p 5 $qq -d $INDIR FM_SNL_09.fds
$QFDS -r -p 5 $qq -d $INDIR FM_SNL_10.fds
$QFDS -r -p 5 $qq -d $INDIR FM_SNL_11.fds
$QFDS -r -p 5 $qq -d $INDIR FM_SNL_12.fds
$QFDS -r -p 5 $qq -d $INDIR FM_SNL_13.fds
$QFDS -r -p 5 $qq -d $INDIR FM_SNL_14.fds
$QFDS -r -p 5 $qq -d $INDIR FM_SNL_15.fds
$QFDS -r -p 5 $qq -d $INDIR FM_SNL_16.fds
$QFDS -r -p 5 $qq -d $INDIR FM_SNL_17.fds
$QFDS -r -p 5 $qq -d $INDIR FM_SNL_21.fds
$QFDS -r -p 5 $qq -d $INDIR FM_SNL_22.fds
