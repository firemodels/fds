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

$QFDS -p 8 -r $qq -d $INDIR Prairie_Grass_Case_21.fds
$QFDS -p 8 -r $qq -d $INDIR Prairie_Grass_Case_22.fds
$QFDS -p 8 -r $qq -d $INDIR Prairie_Grass_Case_23.fds
$QFDS -p 8 -r $qq -d $INDIR Prairie_Grass_Case_24.fds
$QFDS -p 8 -r $qq -d $INDIR Prairie_Grass_Case_42.fds
$QFDS -p 8 -r $qq -d $INDIR Prairie_Grass_Case_46.fds
$QFDS -p 8 -r $qq -d $INDIR Prairie_Grass_Case_55.fds
$QFDS -p 8 -r $qq -d $INDIR Prairie_Grass_Case_57.fds 

echo FDS cases submitted
