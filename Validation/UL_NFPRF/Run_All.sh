#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS  $QUEUE -d $INDIR UL_NFPRF_1_01.fds
$QFDS  $QUEUE -d $INDIR UL_NFPRF_1_02.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_1_03.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_1_04.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_1_05.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_1_06.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_1_07.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_1_08.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_1_09.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_1_10.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_1_11.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_1_12.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_1_13.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_1_14.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_1_15.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_1_16.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_1_17.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_1_18.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_1_19.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_1_20.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_1_21.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_1_22.fds 

$QFDS  $QUEUE -d $INDIR UL_NFPRF_2_01.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_2_02.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_2_03.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_2_04.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_2_05.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_2_06.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_2_07.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_2_08.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_2_09.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_2_10.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_2_11.fds 
$QFDS  $QUEUE -d $INDIR UL_NFPRF_2_12.fds 

$QFDS  $QUEUE -p 16 -d $INDIR UL_NFPRF_P-1.fds
$QFDS  $QUEUE -p 16 -d $INDIR UL_NFPRF_P-2.fds
$QFDS  $QUEUE -p 16 -d $INDIR UL_NFPRF_P-3.fds
$QFDS  $QUEUE -p 16 -d $INDIR UL_NFPRF_P-4.fds
$QFDS  $QUEUE -p 16 -d $INDIR UL_NFPRF_P-5.fds

echo FDS cases submitted
