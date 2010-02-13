#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results McCaffrey_14_kW fire78 &
$RUNFDS Current_Results McCaffrey_22_kW fire78 &
$RUNFDS Current_Results McCaffrey_33_kW fire78 &
$RUNFDS Current_Results McCaffrey_45_kW fire78 &
$RUNFDS Current_Results McCaffrey_57_kW fire78 &

echo FDS cases submitted
