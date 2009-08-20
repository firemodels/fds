@echo off

IF not EXIST placeholder.txt goto dircheck
echo ***error: This script is running in wrong directory.
pause
exit
:dircheck

echo.
echo FDS and Smokeview uninstall script
echo.
echo This script removes path entries for the current FDS/Smokeview installation
echo.
echo Press any key to proceed or CTRL C to abort
pause>NUL

echo.  
echo Proceeding...

echo.
echo removing association between .smv and Smokeview

assoc .smv=
Rem ftype smvDoc=

echo. 
echo Removing FDS and Smokeview shortcuts from the Start menu.
rmdir /q /s "%USERPROFILE%\Start Menu\Programs\FDS5"


cd ..
echo.
echo Removing %CD%\bin from the User Path for: %USERNAME%
call Uninstall\set_path.exe -c "%CD%\bin"

echo.
echo Uninstall Complete.  Press any key to continue.
pause>NUL

