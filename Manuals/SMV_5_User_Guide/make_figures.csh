#!/bin/csh -f
set smokeview=~/bin/smv5_linux
smokeview -help figures/smokeview.help
smokeview -version figures/smokeview.version
smokezip -help figures/smokezip.help
cd ../../Test_cases/Visualization
$smokeview -runscript thouse5
$smokeview -runscript plume5c
$smokeview -runscript smoke_test
$smokeview -runscript smoke_test2
$smokeview -runscript smoke_sensor
$smokeview -runscript sillytexture
$smokeview -runscript plume5a
