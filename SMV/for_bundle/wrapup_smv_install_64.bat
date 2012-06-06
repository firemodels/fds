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
echo Associating the .smv file extension with smokeview_win_64.exe

ftype smvDoc="%CD%\smokeview_win_64.exe" "%%1" >Nul
assoc .smv=smvDoc>Nul

erase "%CD%"\wrapup_smv_install.bat

