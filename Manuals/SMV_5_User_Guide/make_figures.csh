#!/bin/csh -f
echo starting script
set smokeview=~/bin/smv5_linux
set smokezip=~/bin/smokezip_intel
$smokeview -help figures/smokeview.help
$smokeview -version figures/smokeview.version
$smokezip -help figures/smokezip.help
cd ../../Test_cases/Visualization
$smokeview -runscript colorconv
$smokeview -runscript plume5a
$smokeview -runscript plume5b
$smokeview -runscript plume5c
$smokeview -runscript sillytexture
$smokeview -runscript smoke_sensor
$smokeview -runscript smoke_test
$smokeview -runscript smoke_test2
$smokeview -runscript thouse5
$smokeview -runscript script_test
