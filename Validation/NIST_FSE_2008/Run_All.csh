#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results ISOHept19     fire41 &
$RUNFDS Current_Results ISOHept22     fire41 &
$RUNFDS Current_Results ISOHept23     fire44 &
$RUNFDS Current_Results ISOHept24     fire45 &
$RUNFDS Current_Results ISOHept25     fire45 &
$RUNFDS Current_Results ISOHept26     fire46 &
$RUNFDS Current_Results ISOHept27     fire47 &
$RUNFDS Current_Results ISOHept28     fire55 &
$RUNFDS Current_Results ISOHept4      fire51 &
$RUNFDS Current_Results ISOHept5      fire57 &
$RUNFDS Current_Results ISOHept8      fire57 &
$RUNFDS Current_Results ISOHept9      fire58 &
$RUNFDS Current_Results ISOHeptD12    fire59 &
$RUNFDS Current_Results ISOHeptD13    fire59 &
$RUNFDS Current_Results ISONG1        fire62 &
$RUNFDS Current_Results ISONG2        fire62 &
$RUNFDS Current_Results ISONG32       fire62 &
$RUNFDS Current_Results ISONG3        fire65 &
$RUNFDS Current_Results ISONylon10    fire65 &
$RUNFDS Current_Results ISOPP11       fire66 &
$RUNFDS Current_Results ISOPP18       fire66 &
$RUNFDS Current_Results ISOProp15     fire67 &
$RUNFDS Current_Results ISOPropanol30 fire67 &
$RUNFDS Current_Results ISOPropD14    fire68 &
$RUNFDS Current_Results ISOStyrene16  fire68 &
$RUNFDS Current_Results ISOStyrene17  fire70 &
$RUNFDS Current_Results ISOStyrene21  fire70 &
$RUNFDS Current_Results ISOToluene20  fire70 &
$RUNFDS Current_Results ISOToluene29  fire71 &

echo FDS cases submitted
