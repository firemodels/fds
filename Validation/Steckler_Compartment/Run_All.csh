#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results Steckler_010 fire41 &
$RUNFDS Current_Results Steckler_011 fire43 &
$RUNFDS Current_Results Steckler_012 fire44 &
$RUNFDS Current_Results Steckler_013 fire45 &
$RUNFDS Current_Results Steckler_014 fire46 &
$RUNFDS Current_Results Steckler_016 fire47 &
$RUNFDS Current_Results Steckler_017 fire52 &
$RUNFDS Current_Results Steckler_018 fire53 &
$RUNFDS Current_Results Steckler_019 fire54 &
$RUNFDS Current_Results Steckler_020 fire55 &
$RUNFDS Current_Results Steckler_021 fire57 &
$RUNFDS Current_Results Steckler_022 fire58 &
$RUNFDS Current_Results Steckler_023 fire61 &
$RUNFDS Current_Results Steckler_030 fire63 &
$RUNFDS Current_Results Steckler_041 fire63 &
$RUNFDS Current_Results Steckler_114 fire65 &
$RUNFDS Current_Results Steckler_116 fire66 &
$RUNFDS Current_Results Steckler_122 fire66 &
$RUNFDS Current_Results Steckler_144 fire67 &
$RUNFDS Current_Results Steckler_160 fire68 &
$RUNFDS Current_Results Steckler_161 fire62 &
$RUNFDS Current_Results Steckler_162 fire62 &
$RUNFDS Current_Results Steckler_163 fire73 &
$RUNFDS Current_Results Steckler_164 fire73 &
$RUNFDS Current_Results Steckler_165 fire73 &
$RUNFDS Current_Results Steckler_166 fire74 &
$RUNFDS Current_Results Steckler_167 fire74 &
$RUNFDS Current_Results Steckler_210 fire74 &
$RUNFDS Current_Results Steckler_212 fire75 &
$RUNFDS Current_Results Steckler_220 fire76 &
$RUNFDS Current_Results Steckler_221 fire76 &
$RUNFDS Current_Results Steckler_224 fire76 &
$RUNFDS Current_Results Steckler_240 fire77 &
$RUNFDS Current_Results Steckler_242 fire77 &
$RUNFDS Current_Results Steckler_310 fire77 &
$RUNFDS Current_Results Steckler_324 fire78 &
$RUNFDS Current_Results Steckler_410 fire78 &
$RUNFDS Current_Results Steckler_510 fire78 &
$RUNFDS Current_Results Steckler_512 fire77 &
$RUNFDS Current_Results Steckler_513 fire77 &
$RUNFDS Current_Results Steckler_514 fire77 &
$RUNFDS Current_Results Steckler_517 fire78 &
$RUNFDS Current_Results Steckler_520 fire78 &
$RUNFDS Current_Results Steckler_521 fire78 &
$RUNFDS Current_Results Steckler_522 fire78 &
$RUNFDS Current_Results Steckler_524 fire73 &
$RUNFDS Current_Results Steckler_540 fire74 &
$RUNFDS Current_Results Steckler_541 fire75 &
$RUNFDS Current_Results Steckler_542 fire76 &
$RUNFDS Current_Results Steckler_544 fire77 &
$RUNFDS Current_Results Steckler_610 fire78 &
$RUNFDS Current_Results Steckler_612 fire79 &
$RUNFDS Current_Results Steckler_622 fire62 &
$RUNFDS Current_Results Steckler_710 fire62 &
$RUNFDS Current_Results Steckler_810 fire72 &

echo FDS cases submitted
