#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results NIST_NRC_01_v5 fire51 &
$RUNFDS Current_Results NIST_NRC_02_v5 fire51 &
$RUNFDS Current_Results NIST_NRC_03_v5 fire52 &
$RUNFDS Current_Results NIST_NRC_04_v5 fire52 &
$RUNFDS Current_Results NIST_NRC_05_v5 fire53 &
$RUNFDS Current_Results NIST_NRC_07_v5 fire53 &
$RUNFDS Current_Results NIST_NRC_08_v5 fire54 &
$RUNFDS Current_Results NIST_NRC_09_v5 fire54 &
$RUNFDS Current_Results NIST_NRC_10_v5 fire55 &
$RUNFDS Current_Results NIST_NRC_13_v5 fire55 &
$RUNFDS Current_Results NIST_NRC_14_v5 fire56 &
$RUNFDS Current_Results NIST_NRC_15_v5 fire56 &
$RUNFDS Current_Results NIST_NRC_16_v5 fire57 &
$RUNFDS Current_Results NIST_NRC_17_v5 fire57 &
$RUNFDS Current_Results NIST_NRC_18_v5 fire58 &

echo FDS cases submitted
