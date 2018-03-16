@echo off
set CURDIR=%CD%

cd ..\..\Verification\scripts
call Run_FDS_cases.bat -geom
cd %CURDIR%

cd ..\..\Verification\scripts
call Make_GEOM_pictures.bat
cd %CURDIR%
pause

