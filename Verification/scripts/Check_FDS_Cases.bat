@echo off

:: $Date$ 
:: $Revision$
:: $Author$

cd ..
set BASEDIR="%CD%"
cd ..
set SVNROOT="%CD%"
cd %BASEDIR%

set QFDS=call %SVNROOT%\FDS\Utilities\Scripts\checkfds.bat

call FDS_Cases.bat
