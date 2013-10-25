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

$QFDS -r       $qq -d $INDIR backward_facing_step_5.fds
$QFDS -r -p 12 $qq -d $INDIR backward_facing_step_10.fds
#$QFDS -r -p 12 $qq -d $INDIR backward_facing_step_20.fds

echo FDS cases submitted
