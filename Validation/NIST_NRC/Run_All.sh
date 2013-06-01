#!/bin/bash -f

export SVNROOT=`pwd`/../..
export QFDS=/usr/local/bin/qfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results
# qq="-q fire80s"
qq=

# uncomment following line to stop all cases
#export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$QFDS -r $qq -d $INDIR NIST_NRC_01.fds
$QFDS -r $qq -d $INDIR NIST_NRC_02.fds 
$QFDS -r $qq -d $INDIR NIST_NRC_03.fds 
$QFDS -r $qq -d $INDIR NIST_NRC_04.fds 
$QFDS -r $qq -d $INDIR NIST_NRC_05.fds 
$QFDS -r $qq -d $INDIR NIST_NRC_07.fds 
$QFDS -r $qq -d $INDIR NIST_NRC_08.fds 
$QFDS -r $qq -d $INDIR NIST_NRC_09.fds 
$QFDS -r $qq -d $INDIR NIST_NRC_10.fds 
$QFDS -r $qq -d $INDIR NIST_NRC_13.fds 
$QFDS -r $qq -d $INDIR NIST_NRC_14.fds 
$QFDS -r $qq -d $INDIR NIST_NRC_15.fds 
$QFDS -r $qq -d $INDIR NIST_NRC_16.fds 
$QFDS -r $qq -d $INDIR NIST_NRC_17.fds 
$QFDS -r $qq -d $INDIR NIST_NRC_18.fds 

echo FDS cases submitted
