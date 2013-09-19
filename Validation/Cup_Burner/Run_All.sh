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

$QFDS -r -p 4 $qq -d $INDIR Cup_C7H16_Ar.fds
$QFDS -r -p 4 $qq -d $INDIR Cup_C7H16_CO2.fds
$QFDS -r -p 4 $qq -d $INDIR Cup_C7H16_He.fds
$QFDS -r -p 4 $qq -d $INDIR Cup_C7H16_N2.fds
$QFDS -r -p 4 $qq -d $INDIR Cup_CH4_Ar.fds
$QFDS -r -p 4 $qq -d $INDIR Cup_CH4_CO2.fds
$QFDS -r -p 4 $qq -d $INDIR Cup_CH4_He.fds
$QFDS -r -p 4 $qq -d $INDIR Cup_CH4_N2.fds

echo FDS cases submitted
