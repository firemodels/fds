@echo off

echo Creating figures for the Smokeview User's and Verification guides

set vis="%CD%\Visualization"
set wui="%CD%\Wui"
set smvug="%CD%\..\Manuals\SMV_5_User_Guide"
set smvvg="%CD%\..\Manuals\SMV_5_Verification_Guide"

cd %smvug%

erase scriptfigures\*.png
erase scriptfigures\*.help
erase scriptfigures\*.version

smokeview -help > scriptfigures\smokeview.help
smokeview -version > scriptfigures\smokeview.version
smokezip -help > scriptfigures\smokezip.help

cd %smvvg%
erase scriptfigures\*.version
erase scriptfigures\*.png
smokeview -version > scriptfigures\smokeview.version

cd %vis%
smokeview -runscript colorconv
smokeview -runscript devices_elem
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

cd %wui%
smokeview -runscript fire_line

