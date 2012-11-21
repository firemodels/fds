#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR Steckler_010 
$RUNFDS $INDIR Steckler_011 
$RUNFDS $INDIR Steckler_012
$RUNFDS $INDIR Steckler_013
$RUNFDS $INDIR Steckler_014 
$RUNFDS $INDIR Steckler_016 
$RUNFDS $INDIR Steckler_017 
$RUNFDS $INDIR Steckler_018 
$RUNFDS $INDIR Steckler_019 
$RUNFDS $INDIR Steckler_020 
$RUNFDS $INDIR Steckler_021 
$RUNFDS $INDIR Steckler_022 
$RUNFDS $INDIR Steckler_023 
$RUNFDS $INDIR Steckler_030 
$RUNFDS $INDIR Steckler_041 
$RUNFDS $INDIR Steckler_114 
$RUNFDS $INDIR Steckler_116 
$RUNFDS $INDIR Steckler_122 
$RUNFDS $INDIR Steckler_144 
$RUNFDS $INDIR Steckler_160 
$RUNFDS $INDIR Steckler_161 
$RUNFDS $INDIR Steckler_162
$RUNFDS $INDIR Steckler_163 
$RUNFDS $INDIR Steckler_164 
$RUNFDS $INDIR Steckler_165 
$RUNFDS $INDIR Steckler_166 
$RUNFDS $INDIR Steckler_167 
$RUNFDS $INDIR Steckler_210 
$RUNFDS $INDIR Steckler_212 
$RUNFDS $INDIR Steckler_220 
$RUNFDS $INDIR Steckler_221 
$RUNFDS $INDIR Steckler_224 
$RUNFDS $INDIR Steckler_240 
$RUNFDS $INDIR Steckler_242 
$RUNFDS $INDIR Steckler_310 
$RUNFDS $INDIR Steckler_324 
$RUNFDS $INDIR Steckler_410 
$RUNFDS $INDIR Steckler_510 
$RUNFDS $INDIR Steckler_512 
$RUNFDS $INDIR Steckler_513 
$RUNFDS $INDIR Steckler_514 
$RUNFDS $INDIR Steckler_517 
$RUNFDS $INDIR Steckler_520 
$RUNFDS $INDIR Steckler_521 
$RUNFDS $INDIR Steckler_522 
$RUNFDS $INDIR Steckler_524 
$RUNFDS $INDIR Steckler_540 
$RUNFDS $INDIR Steckler_541 
$RUNFDS $INDIR Steckler_542 
$RUNFDS $INDIR Steckler_544 
$RUNFDS $INDIR Steckler_610 
$RUNFDS $INDIR Steckler_612 
$RUNFDS $INDIR Steckler_622 
$RUNFDS $INDIR Steckler_710 
$RUNFDS $INDIR Steckler_810

echo FDS cases submitted
