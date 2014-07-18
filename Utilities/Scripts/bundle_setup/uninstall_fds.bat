@echo off

echo.
echo Uninstalling FDS and Smokeview version 6
NET SESSION >nul 2>&1
IF %ERRORLEVEL% NEQ 0 (
    echo *** Error: This script is running as %username%.  It must run as Administrator.
    echo       Run again, after right clicking on this script and selecting "Run as Adminstrator"
    echo       FDS/Smokeview uninstaller aborted.
    pause
    exit
)
echo.
echo Press any key to proceed or CTRL C to abort
pause>NUL

echo.
echo Removing the association between .smv and Smokeview

assoc .smv=
ftype smvDoc=

echo. 
echo Removing FDS from the Start menu.
rmdir /q /s "%ALLUSERSPROFILE%\Start Menu\Programs\FDS6"

echo.
echo Removing smpd
smpd -remove


