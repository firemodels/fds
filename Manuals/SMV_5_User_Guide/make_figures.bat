@echo off
echo creating some figures for the Smokeview User's guide
erase scriptfigures\*.png
erase scriptfigures\*.help
erase scriptfigures\*.version
smokeview -help scriptfigures\smokeview.help
smokeview -version scriptfigures\smokeview.version
smokezip -help scriptfigures\smokezip.help
cd ..\..\Test_cases\Visualization
smokeview -runscript thouse5
smokeview -runscript plume5c
smokeview -runscript smoke_test
smokeview -runscript smoke_test2
smokeview -runscript smoke_sensor
smokeview -runscript sillytexture
smokeview -runscript plume5a