@echo off
echo creating some figures for the Smokeview User's guide
erase scriptfigures\*.png
cd ..\..\Test_cases\Visualization
smokeview -runscript thouse5
smokeview -runscript plume5c
smokeview -runscript smoke_test
smokeview -runscript smoke_test2