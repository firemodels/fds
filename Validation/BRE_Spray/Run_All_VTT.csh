#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_64/fds5_intel_linux_64
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
setenv BASEDIR `pwd`

cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/Current_Results

$RUNFDS Current_Results BRE_Spray_A_1 
$RUNFDS Current_Results BRE_Spray_A_2
$RUNFDS Current_Results BRE_Spray_A_3 
$RUNFDS Current_Results BRE_Spray_A_4 
$RUNFDS Current_Results BRE_Spray_A_5 
$RUNFDS Current_Results BRE_Spray_A_6 
$RUNFDS Current_Results BRE_Spray_A_7 
$RUNFDS Current_Results BRE_Spray_A_8
$RUNFDS Current_Results BRE_Spray_B_1 
$RUNFDS Current_Results BRE_Spray_B_2
$RUNFDS Current_Results BRE_Spray_B_3 
$RUNFDS Current_Results BRE_Spray_B_4 
$RUNFDS Current_Results BRE_Spray_B_5 
$RUNFDS Current_Results BRE_Spray_B_6 
$RUNFDS Current_Results BRE_Spray_B_7 
$RUNFDS Current_Results BRE_Spray_B_8
$RUNFDS Current_Results BRE_Spray_D_1 
$RUNFDS Current_Results BRE_Spray_D_2
$RUNFDS Current_Results BRE_Spray_D_3 
$RUNFDS Current_Results BRE_Spray_D_4 
$RUNFDS Current_Results BRE_Spray_D_5 
$RUNFDS Current_Results BRE_Spray_D_6 
$RUNFDS Current_Results BRE_Spray_D_7 
$RUNFDS Current_Results BRE_Spray_D_8

echo FDS cases submitted
