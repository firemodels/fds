#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results UL_NFPRF_1_01 fire72 &
$RUNFDS Current_Results UL_NFPRF_1_02 fire72 &
$RUNFDS Current_Results UL_NFPRF_1_03 fire72 &
$RUNFDS Current_Results UL_NFPRF_1_04 fire72 &
$RUNFDS Current_Results UL_NFPRF_1_05 fire73 &
$RUNFDS Current_Results UL_NFPRF_1_06 fire73 &
$RUNFDS Current_Results UL_NFPRF_1_07 fire73 &
$RUNFDS Current_Results UL_NFPRF_1_08 fire73 &
$RUNFDS Current_Results UL_NFPRF_1_09 fire74 &
$RUNFDS Current_Results UL_NFPRF_1_10 fire74 &
$RUNFDS Current_Results UL_NFPRF_1_11 fire74 &
$RUNFDS Current_Results UL_NFPRF_1_12 fire74 &
$RUNFDS Current_Results UL_NFPRF_1_13 fire75 &
$RUNFDS Current_Results UL_NFPRF_1_14 fire75 &
$RUNFDS Current_Results UL_NFPRF_1_15 fire75 &
$RUNFDS Current_Results UL_NFPRF_1_16 fire75 &
$RUNFDS Current_Results UL_NFPRF_1_17 fire76 &
$RUNFDS Current_Results UL_NFPRF_1_18 fire76 &
$RUNFDS Current_Results UL_NFPRF_1_19 fire76 &
$RUNFDS Current_Results UL_NFPRF_1_20 fire76 &
$RUNFDS Current_Results UL_NFPRF_1_21 fire77 &
$RUNFDS Current_Results UL_NFPRF_1_22 fire77 &

echo FDS cases submitted
