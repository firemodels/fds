@echo off

IF not EXIST placeholder.txt goto dircheck
echo ***error: This script is running in the wrong directory.
pause
exit
:dircheck

echo.
echo Wrapping up 64 bit Smokeview update
echo.

echo.
echo Associating the .smv file extension with smokeview.exe

ftype smvDoc="%CD%\smokeview.exe" "%%1" >Nul
assoc .smv=smvDoc>Nul

Rem adding path
echo.
echo Adding %CD% to the system path 
call "%CD%\set_path.exe" -s -m -a "%CD%"


echo Press any key to complete update.
pause>NUL
erase "%CD%|\set_path.exe"
erase "%CD%"\wrapup_smv_install.bat

