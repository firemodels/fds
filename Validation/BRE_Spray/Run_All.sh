#!/bin/bash

export SVNROOT=`pwd`/../..
export QFDS=/usr/local/bin/qfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results
# qq="-q fire80s"
#qq=

# uncomment following line to stop all cases
#export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$QFDS -r $qq -d $INDIR BRE_Spray_A_1.fds
$QFDS -r $qq -d $INDIR BRE_Spray_A_2.fds
$QFDS -r $qq -d $INDIR BRE_Spray_A_3.fds 
$QFDS -r $qq -d $INDIR BRE_Spray_A_4.fds 
$QFDS -r $qq -d $INDIR BRE_Spray_A_5.fds 
$QFDS -r $qq -d $INDIR BRE_Spray_A_6.fds 
$QFDS -r $qq -d $INDIR BRE_Spray_A_7.fds 
$QFDS -r $qq -d $INDIR BRE_Spray_A_8.fds
$QFDS -r $qq -d $INDIR BRE_Spray_B_1.fds 
$QFDS -r $qq -d $INDIR BRE_Spray_B_2.fds
$QFDS -r $qq -d $INDIR BRE_Spray_B_3.fds 
$QFDS -r $qq -d $INDIR BRE_Spray_B_4.fds 
$QFDS -r $qq -d $INDIR BRE_Spray_B_5.fds 
$QFDS -r $qq -d $INDIR BRE_Spray_B_6.fds 
$QFDS -r $qq -d $INDIR BRE_Spray_B_7.fds 
$QFDS -r $qq -d $INDIR BRE_Spray_B_8.fds
$QFDS -r $qq -d $INDIR BRE_Spray_D_1.fds 
$QFDS -r $qq -d $INDIR BRE_Spray_D_2.fds
$QFDS -r $qq -d $INDIR BRE_Spray_D_3.fds 
$QFDS -r $qq -d $INDIR BRE_Spray_D_4.fds 
$QFDS -r $qq -d $INDIR BRE_Spray_D_5.fds 
$QFDS -r $qq -d $INDIR BRE_Spray_D_6.fds 
$QFDS -r $qq -d $INDIR BRE_Spray_D_7.fds 
$QFDS -r $qq -d $INDIR BRE_Spray_D_8.fds
 
echo FDS cases submitted
