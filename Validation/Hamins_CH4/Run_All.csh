#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results Hamins_CH4_01 fire41 &
$RUNFDS Current_Results Hamins_CH4_05 fire42 &
$RUNFDS Current_Results Hamins_CH4_07 fire43 &
$RUNFDS Current_Results Hamins_CH4_19 fire44 &
$RUNFDS Current_Results Hamins_CH4_21 fire45 &
$RUNFDS Current_Results Hamins_CH4_23 fire46 &

echo FDS cases submitted
