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
cd Visualization
$RUNSMV -runscript colorconv
$RUNSMV -runscript devices_elem
$RUNSMV -runscript devices_vistest
$RUNSMV -runscript devices_vistest2
$RUNSMV -runscript plume5a
$RUNSMV -runscript plume5b
$RUNSMV -runscript plume5c
$RUNSMV -runscript plume5c_bounddef
$RUNSMV -runscript sillytexture
$RUNSMV -runscript smoke_sensor
$RUNSMV -runscript smoke_test
$RUNSMV -runscript smoke_test2
$RUNSMV -runscript thouse5
$RUNSMV -runscript script_test
cd ..
cd Wui
$RUNSMV -runscript fire_line


