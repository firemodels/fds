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

$QFDS -r $qq -d $INDIR hu_1p.fds
$QFDS -r $qq -d $INDIR hu_2p.fds
$QFDS -r $qq -d $INDIR hu_3p.fds
$QFDS -r $qq -d $INDIR hu_4p.fds

echo FDS cases submitted
