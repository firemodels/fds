#!/bin/bash -f

export SVNROOT=`pwd`/../..
export QFDS=/usr/local/bin/qfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results
# qq="-q fire80s"
qq=
source ~/.bashrc_fds intel64

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$QFDS -r $qq -d $INDIR Hu_1.fds
$QFDS -r $qq -d $INDIR Hu_2.fds
$QFDS -r $qq -d $INDIR Hu_3.fds
$QFDS -r $qq -d $INDIR Hu_4.fds

echo FDS cases submitted
