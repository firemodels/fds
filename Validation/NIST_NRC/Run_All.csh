#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results NIST_NRC_01 fire51 &
$RUNFDS Current_Results NIST_NRC_02 fire51 &
$RUNFDS Current_Results NIST_NRC_03 fire52 &
$RUNFDS Current_Results NIST_NRC_04 fire52 &
$RUNFDS Current_Results NIST_NRC_05 fire53 &
$RUNFDS Current_Results NIST_NRC_07 fire53 &
$RUNFDS Current_Results NIST_NRC_08 fire54 &
$RUNFDS Current_Results NIST_NRC_09 fire54 &
$RUNFDS Current_Results NIST_NRC_10 fire55 &
$RUNFDS Current_Results NIST_NRC_13 fire55 &
$RUNFDS Current_Results NIST_NRC_14 fire56 &
$RUNFDS Current_Results NIST_NRC_15 fire56 &
$RUNFDS Current_Results NIST_NRC_16 fire57 &
$RUNFDS Current_Results NIST_NRC_17 fire57 &
$RUNFDS Current_Results NIST_NRC_18 fire58 &

echo FDS cases submitted
