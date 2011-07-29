#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results

#uncomment following line to stop all cases
#export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR SE4
$RUNFDS $INDIR SE5
$RUNFDS $INDIR SE6
$RUNFDS $INDIR SE7
$RUNFDS $INDIR SE8
$RUNFDS $INDIR SE9
$RUNFDS $INDIR SE10
$RUNFDS $INDIR SE11
$RUNFDS $INDIR SE12
$RUNFDS $INDIR SE13
$RUNFDS $INDIR SE14
$RUNFDS $INDIR SE15
$RUNFDS $INDIR SE16
$RUNFDS $INDIR SE17
$RUNFDS $INDIR SE18
$RUNFDS $INDIR SE19
$RUNFDS $INDIR SE20
$RUNFDS $INDIR SE21


echo FDS cases submitted
