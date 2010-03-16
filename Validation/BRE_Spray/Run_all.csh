#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`
cp $BASEDIR/FDS_Input_Files/sge-fds-array.sh $BASEDIR/Current_Results

$RUNFDS Current_Results BRE_Spray_1 fire47 &
$RUNFDS Current_Results BRE_Spray_2 fire47 &
$RUNFDS Current_Results BRE_Spray_3 fire47 &
$RUNFDS Current_Results BRE_Spray_4 fire47 &
$RUNFDS Current_Results BRE_Spray_5 fire47 &
$RUNFDS Current_Results BRE_Spray_6 fire47 &
$RUNFDS Current_Results BRE_Spray_7 fire47 &
$RUNFDS Current_Results BRE_Spray_8 fire47 &
$RUNFDS Current_Results BRE_Spray_9 fire47 &
$RUNFDS Current_Results BRE_Spray_10 fire47 &
$RUNFDS Current_Results BRE_Spray_11 fire47 &
$RUNFDS Current_Results BRE_Spray_12 fire47 &
$RUNFDS Current_Results BRE_Spray_13 fire47 &
$RUNFDS Current_Results BRE_Spray_14 fire47 &
$RUNFDS Current_Results BRE_Spray_15 fire47 &
$RUNFDS Current_Results BRE_Spray_16 fire47 &
$RUNFDS Current_Results BRE_Spray_17 fire47 &
$RUNFDS Current_Results BRE_Spray_18 fire47 &
$RUNFDS Current_Results BRE_Spray_19 fire47 &
$RUNFDS Current_Results BRE_Spray_20 fire47 &
$RUNFDS Current_Results BRE_Spray_21 fire47 &
$RUNFDS Current_Results BRE_Spray_22 fire47 &
$RUNFDS Current_Results BRE_Spray_23 fire47 &
$RUNFDS Current_Results BRE_Spray_24 fire47 &

echo FDS cases submitted
