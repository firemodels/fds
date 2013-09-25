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

$QFDS -r $qq -d $INDIR McCaffrey_14_kW.fds 
$QFDS -r $qq -d $INDIR McCaffrey_22_kW.fds 
$QFDS -r $qq -d $INDIR McCaffrey_33_kW.fds 
$QFDS -r $qq -d $INDIR McCaffrey_45_kW.fds 
$QFDS -r $qq -d $INDIR McCaffrey_57_kW.fds 

$QFDS -r $qq -d $INDIR McCaffrey_14_kW_coarse.fds
$QFDS -r $qq -d $INDIR McCaffrey_22_kW_coarse.fds
$QFDS -r $qq -d $INDIR McCaffrey_33_kW_coarse.fds
$QFDS -r $qq -d $INDIR McCaffrey_45_kW_coarse.fds
$QFDS -r $qq -d $INDIR McCaffrey_57_kW_coarse.fds

$QFDS -p 4 -r $qq -d $INDIR McCaffrey_14_kW_fine.fds
$QFDS -p 4 -r $qq -d $INDIR McCaffrey_22_kW_fine.fds
$QFDS -p 4 -r $qq -d $INDIR McCaffrey_33_kW_fine.fds
$QFDS -p 4 -r $qq -d $INDIR McCaffrey_45_kW_fine.fds
$QFDS -p 4 -r $qq -d $INDIR McCaffrey_57_kW_fine.fds

echo FDS cases submitted
