#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results csmag0_32 fire41 &
$RUNFDS Current_Results csmag_32  fire41 &
$RUNFDS Current_Results csmag_64  fire42 &
$RUNFDS Current_Results dsmag_32  fire42 &
$RUNFDS Current_Results dsmag_64  fire43 &
$RUNFDS Current_Results mu0_32    fire43 &

echo FDS cases submitted
