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

$QFDS -r $qq -d $INDIR NIST_Smoke_Alarms_SDC02.fds
$QFDS -r $qq -d $INDIR NIST_Smoke_Alarms_SDC05.fds
$QFDS -r $qq -d $INDIR NIST_Smoke_Alarms_SDC07.fds
$QFDS -r $qq -d $INDIR NIST_Smoke_Alarms_SDC10.fds
$QFDS -r $qq -d $INDIR NIST_Smoke_Alarms_SDC33.fds
$QFDS -r $qq -d $INDIR NIST_Smoke_Alarms_SDC35.fds
$QFDS -r $qq -d $INDIR NIST_Smoke_Alarms_SDC38.fds
$QFDS -r $qq -d $INDIR NIST_Smoke_Alarms_SDC39.fds

echo FDS cases submitted
