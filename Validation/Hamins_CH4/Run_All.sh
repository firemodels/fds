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

$QFDS -r $qq -d $INDIR Hamins_CH4_01.fds
$QFDS -r $qq -d $INDIR Hamins_CH4_05.fds 
$QFDS -r $qq -d $INDIR Hamins_CH4_07.fds 
$QFDS -r $qq -d $INDIR Hamins_CH4_19.fds 
$QFDS -r $qq -d $INDIR Hamins_CH4_21.fds 
$QFDS -r $qq -d $INDIR Hamins_CH4_23.fds 

echo FDS cases submitted
