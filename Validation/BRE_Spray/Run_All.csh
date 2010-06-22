#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`
cp $BASEDIR/FDS_Input_Files/sge-fds.sh $BASEDIR/Current_Results

$RUNFDS Current_Results BRE_Spray_A_1 fire52 &
$RUNFDS Current_Results BRE_Spray_A_2 fire52 &
$RUNFDS Current_Results BRE_Spray_A_3 fire52 &
$RUNFDS Current_Results BRE_Spray_A_4 fire53 &
$RUNFDS Current_Results BRE_Spray_A_5 fire53 &
$RUNFDS Current_Results BRE_Spray_A_6 fire53 &
$RUNFDS Current_Results BRE_Spray_A_7 fire54 &
$RUNFDS Current_Results BRE_Spray_A_8 fire54 &
$RUNFDS Current_Results BRE_Spray_B_1 fire54 &
$RUNFDS Current_Results BRE_Spray_B_2 fire55 &
$RUNFDS Current_Results BRE_Spray_B_3 fire55 &
$RUNFDS Current_Results BRE_Spray_B_4 fire55 &
$RUNFDS Current_Results BRE_Spray_B_5 fire61 &
$RUNFDS Current_Results BRE_Spray_B_6 fire61 &
$RUNFDS Current_Results BRE_Spray_B_7 fire61 &
$RUNFDS Current_Results BRE_Spray_B_8 fire63 &
$RUNFDS Current_Results BRE_Spray_D_1 fire63 &
$RUNFDS Current_Results BRE_Spray_D_2 fire63 &
$RUNFDS Current_Results BRE_Spray_D_3 fire44 &
$RUNFDS Current_Results BRE_Spray_D_4 fire44 &
$RUNFDS Current_Results BRE_Spray_D_5 fire44 &
$RUNFDS Current_Results BRE_Spray_D_6 fire46 &
$RUNFDS Current_Results BRE_Spray_D_7 fire46 &
$RUNFDS Current_Results BRE_Spray_D_8 fire46 &

echo FDS cases submitted
