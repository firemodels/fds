#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results FM_Parallel_Panel_1 fire65 &
$RUNFDS Current_Results FM_Parallel_Panel_2 fire65 &
$RUNFDS Current_Results FM_Parallel_Panel_3 fire66 &
$RUNFDS Current_Results FM_Parallel_Panel_4 fire66 &
$RUNFDS Current_Results FM_Parallel_Panel_5 fire67 &
$RUNFDS Current_Results FM_Parallel_Panel_6 fire67 &

echo FDS cases submitted
