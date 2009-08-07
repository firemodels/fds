#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set RUNSMV=$SVNROOT/SMV_5/Build/INTEL_LINUX_32/smokeview
cd Detectors
$RUNSMV -runscript beam_detector
cd ..
cd Fires
$RUNSMV -runscript room_fire
cd ..
cd Flowfields
$RUNSMV -runscript helium_2d
$RUNSMV -runscript sawtooth
cd ..
cd Miscellaneous
$RUNSMV -runscript pyramid
cd ..
cd NS_Analytical_Solution
$RUNSMV -runscript ns2d_64
cd ..
cd Pressure_Effects
$RUNSMV -runscript pressure_boundary
cd ..
cd Sprinklers_and_Sprays
$RUNSMV -runscript cascade
cd ..


