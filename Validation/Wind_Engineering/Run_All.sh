#/bin/bash

# This script runs the a set of Validation Cases on a linux machine with
# a batch queuing system

export SVNROOT=`pwd`/../..
export QFDS=/usr/local/bin/qfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results
# qq="-q fire80s"
qq=

# uncomment following line to stop all cases
#export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$QFDS -r $qq -p 8 -d $INDIR UWO_test7_case1_180_32.fds
$QFDS -r $qq -p 8 -d $INDIR UWO_test7_case1_180_64.fds
$QFDS -r $qq -p 8 -d $INDIR UWO_test7_case1_270_32.fds
$QFDS -r $qq -p 8 -d $INDIR UWO_test7_case1_270_64.fds
#$QFDS -r $qq -p 8 -d $INDIR UWO_test7_case1_32_spectral.fds
#$QFDS -r $qq -p 8 -d $INDIR UWO_test7_case1_64_spectral.fds

echo FDS cases submitted
