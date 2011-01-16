#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
setenv BASEDIR `pwd`

cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/Current_Results

$RUNFDS Current_Results NRL_HAI_1
$RUNFDS Current_Results NRL_HAI_2
$RUNFDS Current_Results NRL_HAI_3 
$RUNFDS Current_Results NRL_HAI_4
$RUNFDS Current_Results NRL_HAI_5 
$RUNFDS Current_Results NRL_HAI_6 
$RUNFDS Current_Results NRL_HAI_7 
$RUNFDS Current_Results NRL_HAI_8 
$RUNFDS Current_Results NRL_HAI_9 

echo FDS cases submitted
