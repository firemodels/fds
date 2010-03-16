#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`
cp $BASEDIR/FDS_Input_Files/sge-fds.sh $BASEDIR/Current_Results

$RUNFDS Current_Results BRE_Spray_A_1 fire47 &
$RUNFDS Current_Results BRE_Spray_A_2 fire47 &
$RUNFDS Current_Results BRE_Spray_A_3 fire47 &
$RUNFDS Current_Results BRE_Spray_A_4 fire47 &
$RUNFDS Current_Results BRE_Spray_A_5 fire47 &
$RUNFDS Current_Results BRE_Spray_A_6 fire47 &
$RUNFDS Current_Results BRE_Spray_A_7 fire47 &
$RUNFDS Current_Results BRE_Spray_A_8 fire47 &
$RUNFDS Current_Results BRE_Spray_B_1 fire47 &
$RUNFDS Current_Results BRE_Spray_B_2 fire47 &
$RUNFDS Current_Results BRE_Spray_B_3 fire47 &
$RUNFDS Current_Results BRE_Spray_B_4 fire47 &
$RUNFDS Current_Results BRE_Spray_B_5 fire47 &
$RUNFDS Current_Results BRE_Spray_B_6 fire47 &
$RUNFDS Current_Results BRE_Spray_B_7 fire47 &
$RUNFDS Current_Results BRE_Spray_B_8 fire47 &
$RUNFDS Current_Results BRE_Spray_D_1 fire47 &
$RUNFDS Current_Results BRE_Spray_D_2 fire47 &
$RUNFDS Current_Results BRE_Spray_D_3 fire47 &
$RUNFDS Current_Results BRE_Spray_D_4 fire47 &
$RUNFDS Current_Results BRE_Spray_D_5 fire47 &
$RUNFDS Current_Results BRE_Spray_D_6 fire47 &
$RUNFDS Current_Results BRE_Spray_D_7 fire47 &
$RUNFDS Current_Results BRE_Spray_D_8 fire47 &

echo FDS cases submitted
