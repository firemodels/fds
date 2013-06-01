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

$QFDS -r -p 4 $qq -d $INDIR SE4.fds
$QFDS -r -p 4 $qq -d $INDIR SE5.fds
$QFDS -r -p 4 $qq -d $INDIR SE6.fds
$QFDS -r -p 4 $qq -d $INDIR SE7.fds
$QFDS -r -p 4 $qq -d $INDIR SE8.fds
$QFDS -r -p 4 $qq -d $INDIR SE9.fds
$QFDS -r -p 4 $qq -d $INDIR SE10.fds
$QFDS -r -p 4 $qq -d $INDIR SE11.fds
$QFDS -r -p 4 $qq -d $INDIR SE12.fds
$QFDS -r -p 4 $qq -d $INDIR SE13.fds
$QFDS -r -p 4 $qq -d $INDIR SE14.fds
$QFDS -r -p 4 $qq -d $INDIR SE15.fds
$QFDS -r -p 4 $qq -d $INDIR SE16.fds
$QFDS -r -p 4 $qq -d $INDIR SE17.fds
$QFDS -r -p 4 $qq -d $INDIR SE18.fds
$QFDS -r -p 4 $qq -d $INDIR SE19.fds
$QFDS -r -p 4 $qq -d $INDIR SE20.fds
$QFDS -r -p 4 $qq -d $INDIR SE21.fds

echo FDS cases submitted
