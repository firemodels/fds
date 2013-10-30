#!/bin/bash -f

export SVNROOT=`pwd`/../..
export QFDS=/usr/local/bin/qfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results
qq="-q fire80s"
#qq=
source ~/.bashrc_fds intel64

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$QFDS -r $qq -d $INDIR ISOHept19.fds
$QFDS -r $qq -d $INDIR ISOHept22.fds
$QFDS -r $qq -d $INDIR ISOHept23.fds
$QFDS -r $qq -d $INDIR ISOHept24.fds
$QFDS -r $qq -d $INDIR ISOHept25.fds
$QFDS -r $qq -d $INDIR ISOHept26.fds
$QFDS -r $qq -d $INDIR ISOHept27.fds
$QFDS -r $qq -d $INDIR ISOHept28.fds
$QFDS -r $qq -d $INDIR ISOHept4.fds
$QFDS -r $qq -d $INDIR ISOHept5.fds
$QFDS -r $qq -d $INDIR ISOHept8.fds
$QFDS -r $qq -d $INDIR ISOHept9.fds
$QFDS -r $qq -d $INDIR ISOHeptD12.fds
$QFDS -r $qq -d $INDIR ISOHeptD13.fds
$QFDS -r $qq -d $INDIR ISONG1.fds
$QFDS -r $qq -d $INDIR ISONG2.fds
$QFDS -r $qq -d $INDIR ISONG32.fds
$QFDS -r $qq -d $INDIR ISONG3.fds
$QFDS -r $qq -d $INDIR ISONylon10.fds
$QFDS -r $qq -d $INDIR ISOPP11.fds
$QFDS -r $qq -d $INDIR ISOPP18.fds
$QFDS -r $qq -d $INDIR ISOProp15.fds
$QFDS -r $qq -d $INDIR ISOPropanol30.fds
$QFDS -r $qq -d $INDIR ISOPropD14.fds
$QFDS -r $qq -d $INDIR ISOStyrene16.fds
$QFDS -r $qq -d $INDIR ISOStyrene17.fds
$QFDS -r $qq -d $INDIR ISOStyrene21.fds
$QFDS -r $qq -d $INDIR ISOToluene20.fds
$QFDS -r $qq -d $INDIR ISOToluene29.fds

echo FDS cases submitted
