@echo off
set CURDIR=%CD%

cd ..\..\Verification\scripts
call Run_SMV_cases.bat -geom
cd %CURDIR%

cd ..\..\Verification\scripts
call Make_SMV_pictures.bat -geom
cd %CURDIR%
pause

