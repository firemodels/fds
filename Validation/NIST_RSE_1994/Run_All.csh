#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results NIST_RSE_1994_100 fire77 &
$RUNFDS Current_Results NIST_RSE_1994_150 fire77 &
$RUNFDS Current_Results NIST_RSE_1994_200 fire78 &
$RUNFDS Current_Results NIST_RSE_1994_300 fire78 &
$RUNFDS Current_Results NIST_RSE_1994_400 fire78 &
$RUNFDS Current_Results NIST_RSE_1994_500 fire78 &
$RUNFDS Current_Results NIST_RSE_1994_50  fire61 &
$RUNFDS Current_Results NIST_RSE_1994_600 fire61 &
$RUNFDS Current_Results NIST_RSE_1994_75  fire63 &

echo FDS cases submitted
