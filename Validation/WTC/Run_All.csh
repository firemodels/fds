#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results WTC_01 fire63 &
$RUNFDS Current_Results WTC_02 fire63 &
$RUNFDS Current_Results WTC_03 fire64 &
$RUNFDS Current_Results WTC_04 fire64 &
$RUNFDS Current_Results WTC_05 fire78 &
$RUNFDS Current_Results WTC_06 fire78 &

echo FDS cases submitted
