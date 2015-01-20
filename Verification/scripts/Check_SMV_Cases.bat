@echo off

:: $Date$ 
:: $Revision$
:: $Author$

cd ..
set BASEDIR="%CD%"
cd ..
set SVNROOT="%CD%"
cd %BASEDIR%

set QFDS=call %SVNROOT%\Utilities\Scripts\checkfds.bat
set RUNCFAST=call %SVNROOT%\Utilities\Scripts\checkcfast.bat
set RUNTFDS=call %SVNROOT%\Utilities\Scripts\checkfds.bat

echo.
echo Checking FDS cases
echo.

if "%runonlygeom%" == "1" (
  call %SCRIPT_DIR%\SMV_geom_Cases.bat
) else (
  call %SCRIPT_DIR%\SMV_Cases.bat
  call %SCRIPT_DIR%\SMV_geom_Cases.bat
)
