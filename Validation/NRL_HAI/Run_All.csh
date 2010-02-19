#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results NRL_HAI_1 fire44 &
$RUNFDS Current_Results NRL_HAI_2 fire45 &
$RUNFDS Current_Results NRL_HAI_3 fire45 &
$RUNFDS Current_Results NRL_HAI_4 fire46 &
$RUNFDS Current_Results NRL_HAI_5 fire46 &
$RUNFDS Current_Results NRL_HAI_6 fire53 &
$RUNFDS Current_Results NRL_HAI_7 fire53 &
$RUNFDS Current_Results NRL_HAI_8 fire54 &
$RUNFDS Current_Results NRL_HAI_9 fire54 &

echo FDS cases submitted
