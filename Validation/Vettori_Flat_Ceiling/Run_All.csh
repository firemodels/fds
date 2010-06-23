#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`
cp $BASEDIR/FDS_Input_Files/sge-fds.sh $BASEDIR/Current_Results

$RUNFDS Current_Results Vettori_OBSTRUCTED_CORNER_FAST fire68 &
$RUNFDS Current_Results Vettori_OBSTRUCTED_CORNER_MED  fire68 &
$RUNFDS Current_Results Vettori_OBSTRUCTED_CORNER_SLOW fire67 &
$RUNFDS Current_Results Vettori_OBSTRUCTED_OPEN_FAST   fire66 &
$RUNFDS Current_Results Vettori_OBSTRUCTED_OPEN_MED    fire65 &
$RUNFDS Current_Results Vettori_OBSTRUCTED_OPEN_SLOW   fire64 &
$RUNFDS Current_Results Vettori_OBSTRUCTED_WALL_FAST   fire63 &
$RUNFDS Current_Results Vettori_OBSTRUCTED_WALL_MED    fire61 &
$RUNFDS Current_Results Vettori_OBSTRUCTED_WALL_SLOW   fire61 &
$RUNFDS Current_Results Vettori_SMOOTH_CORNER_FAST     fire62 &
$RUNFDS Current_Results Vettori_SMOOTH_CORNER_MED      fire79 &
$RUNFDS Current_Results Vettori_SMOOTH_CORNER_SLOW     fire79 &
$RUNFDS Current_Results Vettori_SMOOTH_OPEN_FAST       fire78 &
$RUNFDS Current_Results Vettori_SMOOTH_OPEN_MED        fire76 &
$RUNFDS Current_Results Vettori_SMOOTH_OPEN_SLOW       fire75 &
$RUNFDS Current_Results Vettori_SMOOTH_WALL_FAST       fire73 &
$RUNFDS Current_Results Vettori_SMOOTH_WALL_MED        fire71 &
$RUNFDS Current_Results Vettori_SMOOTH_WALL_SLOW       fire70 &

echo FDS cases submitted
