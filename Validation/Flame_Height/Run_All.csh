#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results Qs=10000 fire62 &
$RUNFDS Current_Results Qs=1000  fire62 &
$RUNFDS Current_Results Qs=100   fire62 &
$RUNFDS Current_Results Qs=10    fire65 &
$RUNFDS Current_Results Qs=1     fire65 &
$RUNFDS Current_Results Qs=2000  fire65 &
$RUNFDS Current_Results Qs=200   fire66 &
$RUNFDS Current_Results Qs=20    fire66 &
$RUNFDS Current_Results Qs=2     fire67 &
$RUNFDS Current_Results Qs=5000  fire67 &
$RUNFDS Current_Results Qs=500   fire67 &
$RUNFDS Current_Results Qs=50    fire68 &
$RUNFDS Current_Results Qs=5     fire68 &
$RUNFDS Current_Results Qs=p1    fire68 &
$RUNFDS Current_Results Qs=p2    fire62 &
$RUNFDS Current_Results Qs=p5    fire62 &

echo FDS cases submitted
