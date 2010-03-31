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

echo FDS cases submitted
