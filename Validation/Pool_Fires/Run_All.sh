#!/bin/bash

export SVNROOT=`pwd`/../..
export QFDS=/usr/local/bin/qfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results
# qq="-q fire80s"
qq=

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$QFDS -r $qq -d $INDIR acetone_1_m.fds
$QFDS -r $qq -d $INDIR ethanol_1_m.fds
$QFDS -r $qq -d $INDIR methanol_1_m.fds
$QFDS -r $qq -d $INDIR butane_1_m.fds
$QFDS -r $qq -d $INDIR benzene_1_m.fds
$QFDS -r $qq -d $INDIR heptane_1_m.fds


echo FDS cases submitted
