#!/bin/bash

export SVNROOT=`pwd`/../..
export QFDS=/usr/local/bin/qfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results
# qq="-q fire80s"
qq=

# uncomment following line to stop all cases
#export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$QFDS -r $qq -d $INDIR PRISME_LK_1_Lower.fds
$QFDS -r $qq -d $INDIR PRISME_LK_1_Upper.fds
$QFDS -r $qq -d $INDIR PRISME_LK_2_Lower.fds
$QFDS -r $qq -d $INDIR PRISME_LK_2_Upper.fds
$QFDS -r $qq -d $INDIR PRISME_LK_3.fds      
$QFDS -r $qq -d $INDIR PRISME_LK_4.fds      

echo FDS cases submitted
