#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results Hamins_CH4_01 fire79 &
$RUNFDS Current_Results Hamins_CH4_05 fire79 &
$RUNFDS Current_Results Hamins_CH4_07 fire79 &
$RUNFDS Current_Results Hamins_CH4_19 fire79 &
$RUNFDS Current_Results Hamins_CH4_21 fire79 &
$RUNFDS Current_Results Hamins_CH4_23 fire79 &

echo FDS cases submitted
