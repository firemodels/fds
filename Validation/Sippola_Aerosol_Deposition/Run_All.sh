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

$QFDS -r $qq -d $INDIR Sippola_Test_01.fds
$QFDS -r $qq -d $INDIR Sippola_Test_02.fds
$QFDS -r $qq -d $INDIR Sippola_Test_03.fds
$QFDS -r $qq -d $INDIR Sippola_Test_04.fds
$QFDS -r $qq -d $INDIR Sippola_Test_05.fds
$QFDS -r $qq -d $INDIR Sippola_Test_06.fds
$QFDS -r $qq -d $INDIR Sippola_Test_07.fds
$QFDS -r $qq -d $INDIR Sippola_Test_08.fds
$QFDS -r $qq -d $INDIR Sippola_Test_09.fds
$QFDS -r $qq -d $INDIR Sippola_Test_10.fds
$QFDS -r $qq -d $INDIR Sippola_Test_11.fds
$QFDS -r $qq -d $INDIR Sippola_Test_12.fds
$QFDS -r $qq -d $INDIR Sippola_Test_13.fds
$QFDS -r $qq -d $INDIR Sippola_Test_14.fds
$QFDS -r $qq -d $INDIR Sippola_Test_15.fds
$QFDS -r $qq -d $INDIR Sippola_Test_16.fds

echo FDS cases submitted
