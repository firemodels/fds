#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results FM_Parallel_Panel_1 fire41 &
$RUNFDS Current_Results FM_Parallel_Panel_2 fire43 &
$RUNFDS Current_Results FM_Parallel_Panel_3 fire44 &
$RUNFDS Current_Results FM_Parallel_Panel_4 fire45 &
$RUNFDS Current_Results FM_Parallel_Panel_5 fire46 &
$RUNFDS Current_Results FM_Parallel_Panel_6 fire47 &

echo FDS cases submitted
