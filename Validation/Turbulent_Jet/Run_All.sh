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

$QFDS -r $qq -d $INDIR jet_csmag_dx10cm.fds
$QFDS -r $qq -d $INDIR jet_dsmag_dx10cm.fds
$QFDS -r $qq -d $INDIR jet_deardorff_dx10cm.fds
$QFDS -r $qq -d $INDIR jet_vreman_dx10cm.fds

$QFDS -p 49 -r $qq -d $INDIR jet_csmag_dx5cm.fds
$QFDS -p 49 -r $qq -d $INDIR jet_dsmag_dx5cm.fds
$QFDS -p 49 -r $qq -d $INDIR jet_deardorff_dx5cm.fds
$QFDS -p 49 -r $qq -d $INDIR jet_vreman_dx5cm.fds

echo FDS cases submitted
