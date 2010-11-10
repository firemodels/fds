#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results WTC_01 fire71 &
$RUNFDS Current_Results WTC_02 fire71 &
$RUNFDS Current_Results WTC_03 fire72 &
$RUNFDS Current_Results WTC_04 fire72 &
$RUNFDS Current_Results WTC_05 fire79 &
$RUNFDS Current_Results WTC_06 fire79 &

echo FDS cases submitted
