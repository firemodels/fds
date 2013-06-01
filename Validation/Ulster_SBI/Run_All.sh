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

$QFDS -r $qq -d $INDIR Ulster_SBI_30_kW.fds
$QFDS -r $qq -d $INDIR Ulster_SBI_45_kW.fds
$QFDS -r $qq -d $INDIR Ulster_SBI_60_kW.fds 

echo FDS cases submitted
