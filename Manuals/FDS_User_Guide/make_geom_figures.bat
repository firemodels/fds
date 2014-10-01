@echo off
set CURDIR=%CD%

cd ..\..\Verification\scripts
call Run_geom_cases.bat
cd %CURDIR%

cd ..\..\Verification\scripts
call Make_geom_pictures.bat
cd %CURDIR%
pause

