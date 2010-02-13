#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results PRISME_LK_1_Lower fire41 &
$RUNFDS Current_Results PRISME_LK_1_Upper fire41 &
$RUNFDS Current_Results PRISME_LK_2_Lower fire43 &
$RUNFDS Current_Results PRISME_LK_2_Upper fire43 &
$RUNFDS Current_Results PRISME_LK_3       fire44 &
$RUNFDS Current_Results PRISME_LK_4       fire44 &
echo FDS cases submitted
