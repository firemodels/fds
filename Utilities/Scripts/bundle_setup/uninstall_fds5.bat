@echo off

IF not EXIST placeholder.txt goto dircheck
echo ***error: This script is running in wrong directory.
pause
exit
:dircheck

echo.
echo Uninstall FDS and Smokeview
echo.
echo Press any key to proceed or CTRL C to abort
pause>NUL

echo.
echo Removing the association between .smv and Smokeview

assoc .smv=
Rem ftype smvDoc=

echo. 
echo Removing the FDS entry from the Start menu.
rmdir /q /s "%USERPROFILE%\Start Menu\Programs\FDS5"


cd ..
echo.
echo Removing %CD%\bin from the User Path for: %USERNAME%
call Uninstall\set_path.exe -c "%CD%\bin"

echo.
echo Uninstall Complete.  Press any key to continue.
pause>NUL

