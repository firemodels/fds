@echo off

:: $Date$ 
:: $Revision$
:: $Author$

set BASEDIR="%CD%"
cd ..
set SVNROOT="%CD%"
cd %BASEDIR%

set QFDS=call %SVNROOT%\Utilities\Scripts\checkfds.bat

echo.
echo Checking FDS cases
echo.

call FDS_Cases.bat
