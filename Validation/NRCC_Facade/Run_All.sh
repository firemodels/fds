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

$QFDS -r $qq -d $INDIR NRCC_Facade_Win_1_05_MW.fds 
$QFDS -r $qq -d $INDIR NRCC_Facade_Win_1_06_MW.fds 
$QFDS -r $qq -d $INDIR NRCC_Facade_Win_1_08_MW.fds 
$QFDS -r $qq -d $INDIR NRCC_Facade_Win_2_05_MW.fds 
$QFDS -r $qq -d $INDIR NRCC_Facade_Win_2_06_MW.fds 
$QFDS -r $qq -d $INDIR NRCC_Facade_Win_2_08_MW.fds 
$QFDS -r $qq -d $INDIR NRCC_Facade_Win_2_10_MW.fds 
$QFDS -r $qq -d $INDIR NRCC_Facade_Win_3_05_MW.fds 
$QFDS -r $qq -d $INDIR NRCC_Facade_Win_3_06_MW.fds 
$QFDS -r $qq -d $INDIR NRCC_Facade_Win_3_08_MW.fds 
$QFDS -r $qq -d $INDIR NRCC_Facade_Win_3_10_MW.fds 
$QFDS -r $qq -d $INDIR NRCC_Facade_Win_4_05_MW.fds 
$QFDS -r $qq -d $INDIR NRCC_Facade_Win_4_06_MW.fds 
$QFDS -r $qq -d $INDIR NRCC_Facade_Win_4_08_MW.fds 
$QFDS -r $qq -d $INDIR NRCC_Facade_Win_4_10_MW.fds 
$QFDS -r $qq -d $INDIR NRCC_Facade_Win_5_05_MW.fds 
$QFDS -r $qq -d $INDIR NRCC_Facade_Win_5_06_MW.fds 
$QFDS -r $qq -d $INDIR NRCC_Facade_Win_5_08_MW.fds 
$QFDS -r $qq -d $INDIR NRCC_Facade_Win_5_10_MW.fds 

echo FDS cases submitted
