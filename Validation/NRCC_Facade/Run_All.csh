#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results NRCC_Facade_Win_1_05_MW fire41 &
$RUNFDS Current_Results NRCC_Facade_Win_1_06_MW fire41 &
$RUNFDS Current_Results NRCC_Facade_Win_1_08_MW fire43 &
$RUNFDS Current_Results NRCC_Facade_Win_2_05_MW fire43 &
$RUNFDS Current_Results NRCC_Facade_Win_2_06_MW fire45 &
$RUNFDS Current_Results NRCC_Facade_Win_2_08_MW fire45 &
$RUNFDS Current_Results NRCC_Facade_Win_2_10_MW fire46 &
$RUNFDS Current_Results NRCC_Facade_Win_3_05_MW fire46 &
$RUNFDS Current_Results NRCC_Facade_Win_3_06_MW fire47 &
$RUNFDS Current_Results NRCC_Facade_Win_3_08_MW fire47 &
$RUNFDS Current_Results NRCC_Facade_Win_3_10_MW fire56 &
$RUNFDS Current_Results NRCC_Facade_Win_4_05_MW fire56 &
$RUNFDS Current_Results NRCC_Facade_Win_4_06_MW fire57 &
$RUNFDS Current_Results NRCC_Facade_Win_4_08_MW fire57 &
$RUNFDS Current_Results NRCC_Facade_Win_4_10_MW fire58 &
$RUNFDS Current_Results NRCC_Facade_Win_5_05_MW fire58 &
$RUNFDS Current_Results NRCC_Facade_Win_5_06_MW fire59 &
$RUNFDS Current_Results NRCC_Facade_Win_5_08_MW fire59 &
$RUNFDS Current_Results NRCC_Facade_Win_5_10_MW fire61 &

echo FDS cases submitted
