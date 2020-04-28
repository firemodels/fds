#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS -p 64 -n 8 $DEBUG $QUEUE -d $INDIR FM_15cm_Burner_C2H4_15p2_1cm.fds
$QFDS -p 64 -n 8 $DEBUG $QUEUE -d $INDIR FM_15cm_Burner_C2H4_15p2_2cm.fds
$QFDS -p 64 -n 8 $DEBUG $QUEUE -d $INDIR FM_15cm_Burner_C2H4_15p2_5mm.fds
$QFDS -p 64 -n 8 $DEBUG $QUEUE -d $INDIR FM_15cm_Burner_C2H4_16p8_1cm.fds
$QFDS -p 64 -n 8 $DEBUG $QUEUE -d $INDIR FM_15cm_Burner_C2H4_16p8_2cm.fds
$QFDS -p 64 -n 8 $DEBUG $QUEUE -d $INDIR FM_15cm_Burner_C2H4_16p8_5mm.fds
$QFDS -p 64 -n 8 $DEBUG $QUEUE -d $INDIR FM_15cm_Burner_C2H4_19p0_1cm.fds
$QFDS -p 64 -n 8 $DEBUG $QUEUE -d $INDIR FM_15cm_Burner_C2H4_19p0_2cm.fds
$QFDS -p 64 -n 8 $DEBUG $QUEUE -d $INDIR FM_15cm_Burner_C2H4_19p0_5mm.fds
$QFDS -p 64 -n 8 $DEBUG $QUEUE -d $INDIR FM_15cm_Burner_C2H4_20p9_1cm.fds
$QFDS -p 64 -n 8 $DEBUG $QUEUE -d $INDIR FM_15cm_Burner_C2H4_20p9_2cm.fds
$QFDS -p 64 -n 8 $DEBUG $QUEUE -d $INDIR FM_15cm_Burner_C2H4_20p9_5mm.fds

$QFDS -p 125 -n 8 $DEBUG $QUEUE -d $INDIR FM_15cm_Burner_CH4.fds
$QFDS -p 125 -n 8 $DEBUG $QUEUE -d $INDIR FM_15cm_Burner_C2H4.fds
$QFDS -p 125 -n 8 $DEBUG $QUEUE -d $INDIR FM_15cm_Burner_C3H6.fds
$QFDS -p 125 -n 8 $DEBUG $QUEUE -d $INDIR FM_15cm_Burner_C3H8.fds

echo FDS cases submitted
