#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results Ulster_SBI_30_kW fire42 &
$RUNFDS Current_Results Ulster_SBI_45_kW fire43 &
$RUNFDS Current_Results Ulster_SBI_60_kW fire44 &

echo FDS cases submitted
