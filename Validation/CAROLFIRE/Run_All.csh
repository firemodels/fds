#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results CAROLFIRE_PT_01 fire41 &
$RUNFDS Current_Results CAROLFIRE_PT_02 fire41 &
$RUNFDS Current_Results CAROLFIRE_PT_03 fire41 &
$RUNFDS Current_Results CAROLFIRE_PT_04 fire41 &
$RUNFDS Current_Results CAROLFIRE_PT_05 fire41 &
$RUNFDS Current_Results CAROLFIRE_PT_06 fire41 &
$RUNFDS Current_Results CAROLFIRE_PT_07 fire43 &
$RUNFDS Current_Results CAROLFIRE_PT_08 fire43 &
$RUNFDS Current_Results CAROLFIRE_PT_09 fire43 &
$RUNFDS Current_Results CAROLFIRE_PT_10 fire43 &
$RUNFDS Current_Results CAROLFIRE_PT_11 fire43 &
$RUNFDS Current_Results CAROLFIRE_PT_12 fire43 &
$RUNFDS Current_Results CAROLFIRE_PT_13 fire44 &
$RUNFDS Current_Results CAROLFIRE_PT_14 fire44 &
$RUNFDS Current_Results CAROLFIRE_PT_15 fire44 &
$RUNFDS Current_Results CAROLFIRE_PT_16 fire44 &
$RUNFDS Current_Results CAROLFIRE_PT_17 fire45 &
$RUNFDS Current_Results CAROLFIRE_PT_18 fire45 &
$RUNFDS Current_Results CAROLFIRE_PT_19 fire45 &
$RUNFDS Current_Results CAROLFIRE_PT_20 fire45 &
$RUNFDS Current_Results CAROLFIRE_PT_21 fire46 &
$RUNFDS Current_Results CAROLFIRE_PT_22 fire46 &
$RUNFDS Current_Results CAROLFIRE_PT_23 fire46 &
$RUNFDS Current_Results CAROLFIRE_PT_24 fire46 &
$RUNFDS Current_Results CAROLFIRE_PT_25 fire47 &
$RUNFDS Current_Results CAROLFIRE_PT_26 fire47 &
$RUNFDS Current_Results CAROLFIRE_PT_27 fire47 &
$RUNFDS Current_Results CAROLFIRE_PT_28 fire47 &
$RUNFDS Current_Results CAROLFIRE_PT_29 fire41 &
$RUNFDS Current_Results CAROLFIRE_PT_30 fire43 &
$RUNFDS Current_Results CAROLFIRE_PT_31 fire43 &
$RUNFDS Current_Results CAROLFIRE_PT_62 fire44 &
$RUNFDS Current_Results CAROLFIRE_PT_63 fire45 &
$RUNFDS Current_Results CAROLFIRE_PT_64 fire46 &
$RUNFDS Current_Results CAROLFIRE_PT_65 fire47 &

echo FDS cases submitted
