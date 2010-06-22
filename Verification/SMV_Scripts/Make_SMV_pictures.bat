@echo off

echo Creating figures for the Smokeview User's and Verification guides

set vis="%CD%\..\Visualization"
set wui="%CD%\..\Wui"
set smvug="%CD%\..\..\Manuals\SMV_5_User_Guide"
set smvvg="%CD%\..\..\Manuals\SMV_5_Verification_Guide"

cd %smvug%

erase SCRIPT_FIGURES\*.png
erase SCRIPT_FIGURES\*.help
erase SCRIPT_FIGURES\*.version

smokeview -help > SCRIPT_FIGURES\smokeview.help
smokeview -version > SCRIPT_FIGURES\smokeview.version
smokezip -help > SCRIPT_FIGURES\smokezip.help
smokediff -help > SCRIPT_FIGURES\smokediff.help
smokediff -v > SCRIPT_FIGURES\smokediff.version

cd %smvvg%
erase SCRIPT_FIGURES\*.version
erase SCRIPT_FIGURES\*.png
smokeview -version > SCRIPT_FIGURES\smokeview.version

cd %vis%
smokeview -runscript colorbar
smokeview -runscript colorconv
smokeview -runscript objects_elem
smokeview -runscript objects_static
smokeview -runscript objects_dynamic
smokeview -runscript plume5a
smokeview -runscript plume5b
smokeview -runscript plume5c
smokeview -runscript plume5c_bounddef
smokeview -runscript sillytexture
smokeview -runscript smoke_sensor
smokeview -runscript smoke_test
smokeview -runscript smoke_test2
smokeview -runscript thouse5
smokeview -runscript script_test
smokeview -runscript transparency
smokeview -runscript sprinkler_many

cd %wui%
smokeview -runscript fire_line

cd ..\SMV_scripts