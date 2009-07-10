#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results NIST_He_CaseA fire61 &
$RUNFDS Current_Results NIST_He_CaseB fire61 &
$RUNFDS Current_Results NIST_He_CaseC fire62 &
$RUNFDS Current_Results NIST_He_CaseD fire62 &

echo FDS cases submitted
