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

$QFDS -r $qq -d $INDIR Fleury_1t1_100_kW.fds
$QFDS -r $qq -d $INDIR Fleury_1t1_150_kW.fds
$QFDS -r $qq -d $INDIR Fleury_1t1_200_kW.fds
$QFDS -r $qq -d $INDIR Fleury_1t1_250_kW.fds
$QFDS -r $qq -d $INDIR Fleury_1t1_300_kW.fds
$QFDS -r $qq -d $INDIR Fleury_2t1_100_kW.fds
$QFDS -r $qq -d $INDIR Fleury_2t1_150_kW.fds
$QFDS -r $qq -d $INDIR Fleury_2t1_200_kW.fds
$QFDS -r $qq -d $INDIR Fleury_2t1_250_kW.fds
$QFDS -r $qq -d $INDIR Fleury_2t1_300_kW.fds
$QFDS -r $qq -d $INDIR Fleury_3t1_100_kW.fds
$QFDS -r $qq -d $INDIR Fleury_3t1_150_kW.fds
$QFDS -r $qq -d $INDIR Fleury_3t1_200_kW.fds
$QFDS -r $qq -d $INDIR Fleury_3t1_250_kW.fds
$QFDS -r $qq -d $INDIR Fleury_3t1_300_kW.fds

echo FDS cases submitted
