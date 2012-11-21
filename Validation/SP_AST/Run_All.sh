#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR SP_AST_Test_1 
$RUNFDS $INDIR SP_AST_Test_2
$RUNFDS $INDIR SP_AST_Test_3
$RUNFDS $INDIR SP_AST_Diesel_1p1
$RUNFDS $INDIR SP_AST_Diesel_1p9
$RUNFDS $INDIR SP_AST_Heptane_1p1

echo FDS cases submitted
