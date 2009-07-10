#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results Qs=10000 fire41 &
$RUNFDS Current_Results Qs=1000  fire41 &
$RUNFDS Current_Results Qs=100   fire41 &
$RUNFDS Current_Results Qs=10    fire42 &
$RUNFDS Current_Results Qs=1     fire42 &
$RUNFDS Current_Results Qs=2000  fire42 &
$RUNFDS Current_Results Qs=200   fire43 &
$RUNFDS Current_Results Qs=20    fire43 &
$RUNFDS Current_Results Qs=2     fire43 &
$RUNFDS Current_Results Qs=5000  fire44 &
$RUNFDS Current_Results Qs=500   fire44 &
$RUNFDS Current_Results Qs=50    fire44 &
$RUNFDS Current_Results Qs=5     fire45 &
$RUNFDS Current_Results Qs=p1    fire45 &
$RUNFDS Current_Results Qs=p2    fire45 &
$RUNFDS Current_Results Qs=p5    fire45 &

echo FDS cases submitted
