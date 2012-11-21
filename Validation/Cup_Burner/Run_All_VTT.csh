#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_64/fds5_intel_linux_64
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds_sge.csh
setenv BASEDIR `pwd`

cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/Current_Results

$RUNFDS Current_Results Cup_C7H16_Ar
$RUNFDS Current_Results Cup_C7H16_CO2
$RUNFDS Current_Results Cup_C7H16_He
$RUNFDS Current_Results Cup_C7H16_N2
$RUNFDS Current_Results Cup_CH4_Ar
$RUNFDS Current_Results Cup_CH4_CO2
$RUNFDS Current_Results Cup_CH4_He
$RUNFDS Current_Results Cup_CH4_N2

echo FDS cases submitted
