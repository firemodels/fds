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

$QFDS -r $qq -d $INDIR CHRISTIFIRE_S701_tga_N2_v1.fds
$QFDS -r $qq -d $INDIR CHRISTIFIRE_S701_tga_N2_v2.fds
$QFDS -r $qq -d $INDIR CHRISTIFIRE_S701_tga_air_v1.fds
$QFDS -r $qq -d $INDIR CHRISTIFIRE_S701_tga_air_v2.fds
$QFDS -r $qq -d $INDIR CHRISTIFIRE_S701_mcc_v1.fds
$QFDS -r $qq -d $INDIR CHRISTIFIRE_S701_mcc_v2.fds
$QFDS -r $qq -d $INDIR CHRISTIFIRE_I701_tga_N2_v1.fds
$QFDS -r $qq -d $INDIR CHRISTIFIRE_I701_tga_N2_v2.fds
$QFDS -r $qq -d $INDIR CHRISTIFIRE_I701_mcc_v1.fds
$QFDS -r $qq -d $INDIR CHRISTIFIRE_I701_mcc_v2.fds
$QFDS -r $qq -d $INDIR CHRISTIFIRE_C701_cone_25_v1.fds
$QFDS -r $qq -d $INDIR CHRISTIFIRE_C701_cone_25_v2.fds
$QFDS -r $qq -d $INDIR CHRISTIFIRE_C701_cone_50_v1.fds
$QFDS -r $qq -d $INDIR CHRISTIFIRE_C701_cone_50_v2.fds
$QFDS -r $qq -d $INDIR CHRISTIFIRE_C701_cone_75_v1.fds
$QFDS -r $qq -d $INDIR CHRISTIFIRE_C701_cone_75_v2.fds

echo FDS cases submitted
