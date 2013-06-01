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

$QFDS -r $qq -d $INDIR UL_NIST_Vents_Test_1.fds
$QFDS -r $qq -d $INDIR UL_NIST_Vents_Test_2.fds
$QFDS -r $qq -d $INDIR UL_NIST_Vents_Test_3.fds
$QFDS -r $qq -d $INDIR UL_NIST_Vents_Test_4.fds

echo FDS cases submitted
