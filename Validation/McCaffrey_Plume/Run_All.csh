#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results McCaffrey_14_kW fire51 &
$RUNFDS Current_Results McCaffrey_22_kW fire52 &
$RUNFDS Current_Results McCaffrey_33_kW fire53 &
$RUNFDS Current_Results McCaffrey_45_kW fire54 &
$RUNFDS Current_Results McCaffrey_57_kW fire55 &

echo FDS cases submitted
