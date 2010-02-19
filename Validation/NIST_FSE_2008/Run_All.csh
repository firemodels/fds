#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results ISOHept19     fire62 &
$RUNFDS Current_Results ISOHept22     fire62 &
$RUNFDS Current_Results ISOHept23     fire62 &
$RUNFDS Current_Results ISOHept24     fire62 &
$RUNFDS Current_Results ISOHept25     fire66 &
$RUNFDS Current_Results ISOHept26     fire68 &
$RUNFDS Current_Results ISOHept27     fire68 &
$RUNFDS Current_Results ISOHept28     fire51 &
$RUNFDS Current_Results ISOHept4      fire51 &
$RUNFDS Current_Results ISOHept5      fire52 &
$RUNFDS Current_Results ISOHept8      fire52 &
$RUNFDS Current_Results ISOHept9      fire53 &
$RUNFDS Current_Results ISOHeptD12    fire53 &
$RUNFDS Current_Results ISOHeptD13    fire54 &
$RUNFDS Current_Results ISONG1        fire54 &
$RUNFDS Current_Results ISONG2        fire46 &
$RUNFDS Current_Results ISONG32       fire46 &
$RUNFDS Current_Results ISONG3        fire45 &
$RUNFDS Current_Results ISONylon10    fire45 &
$RUNFDS Current_Results ISOPP11       fire44 &
$RUNFDS Current_Results ISOPP18       fire44 &
$RUNFDS Current_Results ISOProp15     fire41 &
$RUNFDS Current_Results ISOPropanol30 fire41 &
$RUNFDS Current_Results ISOPropD14    fire70 &
$RUNFDS Current_Results ISOStyrene16  fire71 &
$RUNFDS Current_Results ISOStyrene17  fire72 &
$RUNFDS Current_Results ISOStyrene21  fire73 &
$RUNFDS Current_Results ISOToluene20  fire74 &
$RUNFDS Current_Results ISOToluene29  fire75 &

echo FDS cases submitted
