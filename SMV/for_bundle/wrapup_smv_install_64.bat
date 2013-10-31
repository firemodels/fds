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

erase "%CD%\set_path.exe"
echo Press any key to complete update
pause>NUL
