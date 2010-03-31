#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results Qs=10000_RI=05 fire41 &
$RUNFDS Current_Results Qs=1000_RI=05  fire41 &
$RUNFDS Current_Results Qs=100_RI=05   fire41 &
$RUNFDS Current_Results Qs=10_RI=05    fire43 &
$RUNFDS Current_Results Qs=1_RI=05     fire43 &
$RUNFDS Current_Results Qs=2000_RI=05  fire43 &
$RUNFDS Current_Results Qs=200_RI=05   fire44 &
$RUNFDS Current_Results Qs=20_RI=05    fire44 &
$RUNFDS Current_Results Qs=2_RI=05     fire44 &
$RUNFDS Current_Results Qs=5000_RI=05  fire45 &
$RUNFDS Current_Results Qs=500_RI=05   fire45 &
$RUNFDS Current_Results Qs=50_RI=05    fire45 &
$RUNFDS Current_Results Qs=5_RI=05     fire46 &
$RUNFDS Current_Results Qs=p1_RI=05    fire46 &
$RUNFDS Current_Results Qs=p2_RI=05    fire47 &
$RUNFDS Current_Results Qs=p5_RI=05    fire47 &

$RUNFDS Current_Results Qs=10000       fire41 &
$RUNFDS Current_Results Qs=1000        fire41 &
$RUNFDS Current_Results Qs=100         fire41 &
$RUNFDS Current_Results Qs=10          fire43 &
$RUNFDS Current_Results Qs=1           fire43 &
$RUNFDS Current_Results Qs=2000        fire43 &
$RUNFDS Current_Results Qs=200         fire44 &
$RUNFDS Current_Results Qs=20          fire44 &
$RUNFDS Current_Results Qs=2           fire44 &
$RUNFDS Current_Results Qs=5000        fire45 &
$RUNFDS Current_Results Qs=500         fire45 &
$RUNFDS Current_Results Qs=50          fire45 &
$RUNFDS Current_Results Qs=5           fire46 &
$RUNFDS Current_Results Qs=p1          fire46 &
$RUNFDS Current_Results Qs=p2          fire47 &
$RUNFDS Current_Results Qs=p5          fire47 &

$RUNFDS Current_Results Qs=10000_RI=20 fire51 &
$RUNFDS Current_Results Qs=1000_RI=20  fire51 &
$RUNFDS Current_Results Qs=100_RI=20   fire52 &
$RUNFDS Current_Results Qs=10_RI=20    fire52 &
$RUNFDS Current_Results Qs=1_RI=20     fire53 &
$RUNFDS Current_Results Qs=2000_RI=20  fire53 &
$RUNFDS Current_Results Qs=200_RI=20   fire54 &
$RUNFDS Current_Results Qs=20_RI=20    fire54 &
$RUNFDS Current_Results Qs=2_RI=20     fire55 &
$RUNFDS Current_Results Qs=5000_RI=20  fire55 &
$RUNFDS Current_Results Qs=500_RI=20   fire56 &
$RUNFDS Current_Results Qs=50_RI=20    fire56 &
$RUNFDS Current_Results Qs=5_RI=20     fire57 &
$RUNFDS Current_Results Qs=p1_RI=20    fire57 &
$RUNFDS Current_Results Qs=p2_RI=20    fire58 &
$RUNFDS Current_Results Qs=p5_RI=20    fire58 &

echo FDS cases submitted
