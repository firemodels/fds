#!/bin/csh -f

#  Note: this script can only be run at the console of a Linux or OSX computer

#  To use this script:
#  build smokeview first by cd'ing to ../SMV/Build/INTEL_OSX_32 
#  or directory for whatever version of smokeview you wish to build .
#  Then type:  ./make_smv.csh
#  Look in the makefile to see how smokeview name is spelled and 
#  edit line below appropriately .


#set smokeview=../SMV/bin/smv5_linux_32
#set smokeview=../SMV/bin/smv5_linux_64
set smokeview=../SMV/bin/smv5_osx_32
#set smokeview=../SMV/bin/smv5_osx_64

echo
echo erasing User Guide scripted figures generated previously
rm ../Manuals/FDS_User_Guide/SCRIPT_FIGURES/*.png

echo
echo erasing Verification scripted figures generated previously
erase ../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/*.png

cd Controls
$smokeview -runscript activate_vents

cd ../Detectors
$smokeview -runscript beam_detector

cd ../Fires
$smokeview -runscript room_fire

cd ../Flowfields
$smokeview -runscript helium_2d
$smokeview -runscript sawtooth
$smokeview -runscript symmetry_test

cd ../HVAC
$smokeview -runscript HVAC_mass_conservation
$smokeview -runscript HVAC_energy_pressure

cd ../Miscellaneous
$smokeview -runscript pyramid

cd ../NS_Analytical_Solution
$smokeview -runscript ns2d_64

cd ../Pressure_Effects
$smokeview -runscript pressure_boundary
$smokeview -runscript leak_test_2

cd ../Scalar_Analytical_Solution
$smokeview -runscript move_slug
$smokeview -runscript move_slug_fl1

cd ../Sprinklers_and_Sprays
$smokeview -runscript cascade

cd ../Visualization
$smokeview -runscript objects_static
$smokeview -runscript objects_dynamic


# generate Smokeview figures
#cd ../scripts
#./Make_SMV_Pictures


