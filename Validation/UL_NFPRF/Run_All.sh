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

$QFDS -r $qq -d $INDIR UL_NFPRF_1_01.fds
$QFDS -r $qq -d $INDIR UL_NFPRF_1_02.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_1_03.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_1_04.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_1_05.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_1_06.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_1_07.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_1_08.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_1_09.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_1_10.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_1_11.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_1_12.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_1_13.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_1_14.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_1_15.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_1_16.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_1_17.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_1_18.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_1_19.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_1_20.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_1_21.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_1_22.fds 

$QFDS -r $qq -d $INDIR UL_NFPRF_2_01.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_2_02.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_2_03.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_2_04.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_2_05.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_2_06.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_2_07.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_2_08.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_2_09.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_2_10.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_2_11.fds 
$QFDS -r $qq -d $INDIR UL_NFPRF_2_12.fds 

echo FDS cases submitted
