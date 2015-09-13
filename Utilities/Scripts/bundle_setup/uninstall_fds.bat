@echo off

NET SESSION >nul 2>&1
IF %ERRORLEVEL% NEQ 0 (
    echo *** Error: This script is running as %username%.  It must run as Administrator.
    echo       Run again, after right clicking on this script and selecting "Run as Adminstrator"
    echo       FDS/Smokeview uninstaller aborted.
    pause
    exit
)

call :is_cfast_installed
if %cfastinstalled% == 1 goto skip1
  echo.
  echo *** Removing the association between .smv and Smokeview
  assoc .smv=
  ftype smvDoc=
:skip1

echo. 
echo *** Removing FDS from the Start menu.
rmdir /q /s "%ALLUSERSPROFILE%\Start Menu\Programs\FDS6"

echo.
echo *** Removing hydra_service
hydra_service -remove >Nul 2>Nul

echo.
echo *** Removing smpd
smpd -remove >Nul 2>Nul
