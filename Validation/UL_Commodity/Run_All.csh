#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results Group_A_2x2x2 fire54 &
$RUNFDS Current_Results Group_A_2x2x3 fire72 &
$RUNFDS Current_Results Group_A_2x2x4 fire72 &

echo FDS cases submitted
