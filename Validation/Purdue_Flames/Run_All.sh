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

$QFDS -r $qq -d $INDIR 7p1_cm_methane_4mm.fds
#$QFDS -r $qq -d $INDIR 7p1_cm_methane_2mm.fds
$QFDS -r -p 16 $qq -d $INDIR 7p1_cm_methane_2mm_16mesh.fds
#$QFDS -r -p 128 $qq -d $INDIR 7p1_cm_methane_1mm_128mesh_les.fds
#$QFDS -r -p 128 $qq -d $INDIR 7p1_cm_methane_1mm_128mesh_dns.fds

echo FDS cases submitted
