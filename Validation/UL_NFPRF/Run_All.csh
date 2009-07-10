#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results UL_NFPRF_1_01 fire71 &
$RUNFDS Current_Results UL_NFPRF_1_02 fire71 &
$RUNFDS Current_Results UL_NFPRF_1_03 fire72 &
$RUNFDS Current_Results UL_NFPRF_1_04 fire72 &
$RUNFDS Current_Results UL_NFPRF_1_05 fire73 &
$RUNFDS Current_Results UL_NFPRF_1_06 fire73 &
$RUNFDS Current_Results UL_NFPRF_1_07 fire74 &
$RUNFDS Current_Results UL_NFPRF_1_08 fire74 &
$RUNFDS Current_Results UL_NFPRF_1_09 fire75 &
$RUNFDS Current_Results UL_NFPRF_1_10 fire75 &
$RUNFDS Current_Results UL_NFPRF_1_11 fire76 &
$RUNFDS Current_Results UL_NFPRF_1_12 fire76 &
$RUNFDS Current_Results UL_NFPRF_1_13 fire77 &
$RUNFDS Current_Results UL_NFPRF_1_14 fire77 &
$RUNFDS Current_Results UL_NFPRF_1_15 fire78 &
$RUNFDS Current_Results UL_NFPRF_1_16 fire78 &
$RUNFDS Current_Results UL_NFPRF_1_17 fire79 &
$RUNFDS Current_Results UL_NFPRF_1_18 fire79 &
$RUNFDS Current_Results UL_NFPRF_1_19 fire79 &
$RUNFDS Current_Results UL_NFPRF_1_20 fire73 &
$RUNFDS Current_Results UL_NFPRF_1_21 fire74 &
$RUNFDS Current_Results UL_NFPRF_1_22 fire75 &

echo FDS cases submitted
