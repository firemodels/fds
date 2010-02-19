#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results NIST_He_CaseA fire65 &
$RUNFDS Current_Results NIST_He_CaseB fire65 &
$RUNFDS Current_Results NIST_He_CaseC fire67 &
$RUNFDS Current_Results NIST_He_CaseD fire67 &

echo FDS cases submitted
