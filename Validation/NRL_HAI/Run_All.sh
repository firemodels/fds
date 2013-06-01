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

$QFDS -r $qq -d $INDIR NRL_HAI_1.fds
$QFDS -r $qq -d $INDIR NRL_HAI_2.fds
$QFDS -r $qq -d $INDIR NRL_HAI_3.fds 
$QFDS -r $qq -d $INDIR NRL_HAI_4.fds
$QFDS -r $qq -d $INDIR NRL_HAI_5.fds 
$QFDS -r $qq -d $INDIR NRL_HAI_6.fds 
$QFDS -r $qq -d $INDIR NRL_HAI_7.fds 
$QFDS -r $qq -d $INDIR NRL_HAI_8.fds 
$QFDS -r $qq -d $INDIR NRL_HAI_9.fds 

echo FDS cases submitted
