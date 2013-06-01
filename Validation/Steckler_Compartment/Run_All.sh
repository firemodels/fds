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

$QFDS -r $qq -d $INDIR Steckler_010.fds
$QFDS -r $qq -d $INDIR Steckler_011.fds
$QFDS -r $qq -d $INDIR Steckler_012.fds
$QFDS -r $qq -d $INDIR Steckler_013.fds
$QFDS -r $qq -d $INDIR Steckler_014.fds
$QFDS -r $qq -d $INDIR Steckler_016.fds
$QFDS -r $qq -d $INDIR Steckler_017.fds
$QFDS -r $qq -d $INDIR Steckler_018.fds
$QFDS -r $qq -d $INDIR Steckler_019.fds
$QFDS -r $qq -d $INDIR Steckler_020.fds
$QFDS -r $qq -d $INDIR Steckler_021.fds
$QFDS -r $qq -d $INDIR Steckler_022.fds
$QFDS -r $qq -d $INDIR Steckler_023.fds
$QFDS -r $qq -d $INDIR Steckler_030.fds
$QFDS -r $qq -d $INDIR Steckler_041.fds
$QFDS -r $qq -d $INDIR Steckler_114.fds
$QFDS -r $qq -d $INDIR Steckler_116.fds
$QFDS -r $qq -d $INDIR Steckler_122.fds
$QFDS -r $qq -d $INDIR Steckler_144.fds
$QFDS -r $qq -d $INDIR Steckler_160.fds
$QFDS -r $qq -d $INDIR Steckler_161.fds
$QFDS -r $qq -d $INDIR Steckler_162.fds
$QFDS -r $qq -d $INDIR Steckler_163.fds
$QFDS -r $qq -d $INDIR Steckler_164.fds
$QFDS -r $qq -d $INDIR Steckler_165.fds
$QFDS -r $qq -d $INDIR Steckler_166.fds
$QFDS -r $qq -d $INDIR Steckler_167.fds
$QFDS -r $qq -d $INDIR Steckler_210.fds
$QFDS -r $qq -d $INDIR Steckler_212.fds
$QFDS -r $qq -d $INDIR Steckler_220.fds
$QFDS -r $qq -d $INDIR Steckler_221.fds
$QFDS -r $qq -d $INDIR Steckler_224.fds
$QFDS -r $qq -d $INDIR Steckler_240.fds
$QFDS -r $qq -d $INDIR Steckler_242.fds
$QFDS -r $qq -d $INDIR Steckler_310.fds
$QFDS -r $qq -d $INDIR Steckler_324.fds
$QFDS -r $qq -d $INDIR Steckler_410.fds
$QFDS -r $qq -d $INDIR Steckler_510.fds
$QFDS -r $qq -d $INDIR Steckler_512.fds
$QFDS -r $qq -d $INDIR Steckler_513.fds
$QFDS -r $qq -d $INDIR Steckler_514.fds
$QFDS -r $qq -d $INDIR Steckler_517.fds
$QFDS -r $qq -d $INDIR Steckler_520.fds
$QFDS -r $qq -d $INDIR Steckler_521.fds
$QFDS -r $qq -d $INDIR Steckler_522.fds
$QFDS -r $qq -d $INDIR Steckler_524.fds
$QFDS -r $qq -d $INDIR Steckler_540.fds
$QFDS -r $qq -d $INDIR Steckler_541.fds
$QFDS -r $qq -d $INDIR Steckler_542.fds
$QFDS -r $qq -d $INDIR Steckler_544.fds
$QFDS -r $qq -d $INDIR Steckler_610.fds
$QFDS -r $qq -d $INDIR Steckler_612.fds
$QFDS -r $qq -d $INDIR Steckler_622.fds
$QFDS -r $qq -d $INDIR Steckler_710.fds
$QFDS -r $qq -d $INDIR Steckler_810.fds

echo FDS cases submitted
