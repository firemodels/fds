#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results WTC_01_v5 fire61 &
$RUNFDS Current_Results WTC_02_v5 fire61 &
$RUNFDS Current_Results WTC_03_v5 fire63 &
$RUNFDS Current_Results WTC_04_v5 fire63 &
$RUNFDS Current_Results WTC_05_v5 fire64 &
$RUNFDS Current_Results WTC_06_v5 fire64 &

echo FDS cases submitted
