@echo off

set SCRIPT_DIR=%CD%
cd ..
set BASEDIR="%CD%"
cd ..
set SVNROOT="%CD%"
cd %BASEDIR%

set QFDS=call %SVNROOT%\Utilities\Scripts\checkfds.bat
set RUNCFAST=call %SVNROOT%\Utilities\Scripts\checkcfast.bat
set RUNTFDS=call %SVNROOT%\Utilities\Scripts\checkfds.bat

call %SCRIPT_DIR%\SMV_Cases.bat
call %SCRIPT_DIR%\GEOM_Cases.bat
call %SCRIPT_DIR%\WUI_Cases.bat
