#!/bin/bash -f

export SVNROOT=`pwd`/../..
export QFDS=/usr/local/bin/qfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results
# qq="-q fire80s"
qq=

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$QFDS -r $qq -d $INDIR Bryant_034_kW.fds
$QFDS -r $qq -d $INDIR Bryant_065_kW.fds
$QFDS -r $qq -d $INDIR Bryant_096_kW.fds
$QFDS -r $qq -d $INDIR Bryant_128_kW.fds
$QFDS -r $qq -d $INDIR Bryant_160_kW.fds
$QFDS -r $qq -d $INDIR Bryant_320_kW.fds
$QFDS -r $qq -d $INDIR Bryant_511_kW.fds
 
echo FDS cases submitted
