#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results ISOHept19     fire61 &
$RUNFDS Current_Results ISOHept22     fire61 &
$RUNFDS Current_Results ISOHept23     fire62 &
$RUNFDS Current_Results ISOHept24     fire62 &
$RUNFDS Current_Results ISOHept25     fire63 &
$RUNFDS Current_Results ISOHept26     fire63 &
$RUNFDS Current_Results ISOHept27     fire64 &
$RUNFDS Current_Results ISOHept28     fire64 &
$RUNFDS Current_Results ISOHept4      fire65 &
$RUNFDS Current_Results ISOHept5      fire65 &
$RUNFDS Current_Results ISOHept8      fire66 &
$RUNFDS Current_Results ISOHept9      fire66 &
$RUNFDS Current_Results ISOHeptD12    fire67 &
$RUNFDS Current_Results ISOHeptD13    fire67 &
$RUNFDS Current_Results ISONG1        fire68 &
$RUNFDS Current_Results ISONG2        fire68 &
$RUNFDS Current_Results ISONG32       fire70 &
$RUNFDS Current_Results ISONG3        fire70 &
$RUNFDS Current_Results ISONylon10    fire71 &
$RUNFDS Current_Results ISOPP11       fire71 &
$RUNFDS Current_Results ISOPP18       fire73 &
$RUNFDS Current_Results ISOProp15     fire73 &
$RUNFDS Current_Results ISOPropanol30 fire74 &
$RUNFDS Current_Results ISOPropD14    fire74 &
$RUNFDS Current_Results ISOStyrene16  fire75 &
$RUNFDS Current_Results ISOStyrene17  fire75 &
$RUNFDS Current_Results ISOStyrene21  fire76 &
$RUNFDS Current_Results ISOToluene20  fire76 &
$RUNFDS Current_Results ISOToluene29  fire77 &

echo FDS cases submitted
