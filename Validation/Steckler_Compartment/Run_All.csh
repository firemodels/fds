#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

$RUNFDS Current_Results Steckler_010 fire70 &
$RUNFDS Current_Results Steckler_011 fire70 &
$RUNFDS Current_Results Steckler_012 fire70 &
$RUNFDS Current_Results Steckler_013 fire70 &
$RUNFDS Current_Results Steckler_014 fire70 &
$RUNFDS Current_Results Steckler_016 fire70 &
$RUNFDS Current_Results Steckler_017 fire70 &
$RUNFDS Current_Results Steckler_018 fire70 &
$RUNFDS Current_Results Steckler_019 fire71 &
$RUNFDS Current_Results Steckler_020 fire71 &
$RUNFDS Current_Results Steckler_021 fire71 &
$RUNFDS Current_Results Steckler_022 fire71 &
$RUNFDS Current_Results Steckler_023 fire71 &
$RUNFDS Current_Results Steckler_030 fire71 &
$RUNFDS Current_Results Steckler_041 fire71 &
$RUNFDS Current_Results Steckler_114 fire71 &
$RUNFDS Current_Results Steckler_116 fire73 &
$RUNFDS Current_Results Steckler_122 fire73 &
$RUNFDS Current_Results Steckler_144 fire73 &
$RUNFDS Current_Results Steckler_160 fire73 &
$RUNFDS Current_Results Steckler_161 fire73 &
$RUNFDS Current_Results Steckler_162 fire73 &
$RUNFDS Current_Results Steckler_163 fire73 &
$RUNFDS Current_Results Steckler_164 fire73 &
$RUNFDS Current_Results Steckler_165 fire74 &
$RUNFDS Current_Results Steckler_166 fire74 &
$RUNFDS Current_Results Steckler_167 fire74 &
$RUNFDS Current_Results Steckler_210 fire74 &
$RUNFDS Current_Results Steckler_212 fire74 &
$RUNFDS Current_Results Steckler_220 fire74 &
$RUNFDS Current_Results Steckler_221 fire74 &
$RUNFDS Current_Results Steckler_224 fire74 &
$RUNFDS Current_Results Steckler_240 fire75 &
$RUNFDS Current_Results Steckler_242 fire75 &
$RUNFDS Current_Results Steckler_310 fire75 &
$RUNFDS Current_Results Steckler_324 fire75 &
$RUNFDS Current_Results Steckler_410 fire75 &
$RUNFDS Current_Results Steckler_510 fire75 &
$RUNFDS Current_Results Steckler_512 fire75 &
$RUNFDS Current_Results Steckler_513 fire75 &
$RUNFDS Current_Results Steckler_514 fire76 &
$RUNFDS Current_Results Steckler_517 fire76 &
$RUNFDS Current_Results Steckler_520 fire76 &
$RUNFDS Current_Results Steckler_521 fire76 &
$RUNFDS Current_Results Steckler_522 fire76 &
$RUNFDS Current_Results Steckler_524 fire76 &
$RUNFDS Current_Results Steckler_540 fire76 &
$RUNFDS Current_Results Steckler_541 fire76 &
$RUNFDS Current_Results Steckler_542 fire77 &
$RUNFDS Current_Results Steckler_544 fire77 &
$RUNFDS Current_Results Steckler_610 fire77 &
$RUNFDS Current_Results Steckler_612 fire77 &
$RUNFDS Current_Results Steckler_622 fire77 &
$RUNFDS Current_Results Steckler_710 fire77 &
$RUNFDS Current_Results Steckler_810 fire77 &

echo FDS cases submitted
