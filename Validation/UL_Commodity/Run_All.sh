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

$QFDS -r $qq -d $INDIR Group_A_2x2x2.fds
$QFDS -r $qq -d $INDIR Group_A_2x2x3.fds
$QFDS -r $qq -d $INDIR Group_A_2x2x4.fds
$QFDS -r $qq -d $INDIR Group_A_FM_RDD_p21.fds
$QFDS -r $qq  -d $INDIR Group_A_FM_RDD_p31.fds
$QFDS -r $qq  -d $INDIR Group_A_FM_RDD_p39.fds

echo FDS cases submitted
