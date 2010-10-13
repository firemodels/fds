#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results UL_NFPRF_1_01 fire61 &
$RUNFDS Current_Results UL_NFPRF_1_02 fire61 &
$RUNFDS Current_Results UL_NFPRF_1_03 fire62 &
$RUNFDS Current_Results UL_NFPRF_1_04 fire62 &
$RUNFDS Current_Results UL_NFPRF_1_05 fire63 &
$RUNFDS Current_Results UL_NFPRF_1_06 fire63 &
$RUNFDS Current_Results UL_NFPRF_1_07 fire64 &
$RUNFDS Current_Results UL_NFPRF_1_08 fire64 &
$RUNFDS Current_Results UL_NFPRF_1_09 fire65 &
$RUNFDS Current_Results UL_NFPRF_1_10 fire65 &
$RUNFDS Current_Results UL_NFPRF_1_11 fire66 &
$RUNFDS Current_Results UL_NFPRF_1_12 fire66 &
$RUNFDS Current_Results UL_NFPRF_1_13 fire67 &
$RUNFDS Current_Results UL_NFPRF_1_14 fire67 &
$RUNFDS Current_Results UL_NFPRF_1_15 fire68 &
$RUNFDS Current_Results UL_NFPRF_1_16 fire68 &
$RUNFDS Current_Results UL_NFPRF_1_17 fire75 &
$RUNFDS Current_Results UL_NFPRF_1_18 fire75 &
$RUNFDS Current_Results UL_NFPRF_1_19 fire75 &
$RUNFDS Current_Results UL_NFPRF_1_20 fire76 &
$RUNFDS Current_Results UL_NFPRF_1_21 fire76 &
$RUNFDS Current_Results UL_NFPRF_1_22 fire76 &

$RUNFDS Current_Results UL_NFPRF_2_01 fire70 &
$RUNFDS Current_Results UL_NFPRF_2_02 fire70 &
$RUNFDS Current_Results UL_NFPRF_2_03 fire70 &
$RUNFDS Current_Results UL_NFPRF_2_04 fire71 &
$RUNFDS Current_Results UL_NFPRF_2_05 fire71 &
$RUNFDS Current_Results UL_NFPRF_2_06 fire71 &
$RUNFDS Current_Results UL_NFPRF_2_07 fire72 &
$RUNFDS Current_Results UL_NFPRF_2_08 fire72 &
$RUNFDS Current_Results UL_NFPRF_2_09 fire72 &
$RUNFDS Current_Results UL_NFPRF_2_10 fire73 &
$RUNFDS Current_Results UL_NFPRF_2_11 fire73 &
$RUNFDS Current_Results UL_NFPRF_2_12 fire73 &

echo FDS cases submitted
