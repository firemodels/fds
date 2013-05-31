#!/bin/bash -f

export SVNROOT=`pwd`/../..
export QFDS=/usr/local/bin/qfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results
# qq="-q fire80s"
qq=

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$QFDS -r $qq -d $INDIR CAROLFIRE_PT_01.fds
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_02.fds
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_03.fds
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_04.fds
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_05.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_06.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_07.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_08.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_09.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_10.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_11.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_12.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_13.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_14.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_15.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_16.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_17.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_18.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_19.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_20.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_21.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_22.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_23.fds
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_24.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_25.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_26.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_27.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_28.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_29.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_30.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_31.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_62.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_63.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_64.fds 
$QFDS -r $qq -d $INDIR CAROLFIRE_PT_65.fds 

echo FDS cases submitted
