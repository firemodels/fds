#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results Hamins_CH4_01 fire61 &
$RUNFDS Current_Results Hamins_CH4_05 fire61 &
$RUNFDS Current_Results Hamins_CH4_07 fire62 &
$RUNFDS Current_Results Hamins_CH4_19 fire62 &
$RUNFDS Current_Results Hamins_CH4_21 fire63 &
$RUNFDS Current_Results Hamins_CH4_23 fire63 &

echo FDS cases submitted
