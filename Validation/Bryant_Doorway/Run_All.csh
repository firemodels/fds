#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results Bryant_034_kW fire51 &
$RUNFDS Current_Results Bryant_065_kW fire54 &
$RUNFDS Current_Results Bryant_096_kW fire55 &
$RUNFDS Current_Results Bryant_128_kW fire57 &
$RUNFDS Current_Results Bryant_160_kW fire58 &
$RUNFDS Current_Results Bryant_320_kW fire59 &
$RUNFDS Current_Results Bryant_511_kW fire59 &
 
echo FDS cases submitted
