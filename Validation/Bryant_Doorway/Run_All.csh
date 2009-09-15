#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results Bryant_034_kW fire65 &
$RUNFDS Current_Results Bryant_065_kW fire65 &
$RUNFDS Current_Results Bryant_096_kW fire66 &
$RUNFDS Current_Results Bryant_128_kW fire66 &
$RUNFDS Current_Results Bryant_160_kW fire67 &
$RUNFDS Current_Results Bryant_320_kW fire67 &
$RUNFDS Current_Results Bryant_511_kW fire68 &
 
echo FDS cases submitted
