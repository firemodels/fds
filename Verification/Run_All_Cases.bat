@echo off
set svn_drive=d:

set BASEDIR=%CD%
set SVNROOT=%BASEDIR%\..\

set FDS=%SVNROOT%\FDS_Compilation\intel_win_32\fds5_win_32

set RUNFDS=call %SVNROOT%\Utilities\Scripts\runfds_win32.bat

echo You are about to run the Verification Test Suite.
pause > Nul

echo FDS cases submitted

call FDS_Cases.bat

cd %BASEDIR%

Rem create a text file containing the FDS5 version used to run these tests.
Rem This file is included in the smokeview user's guide

set smvug="%CD%\..\Manuals\SMV_5_User_Guide\"
echo | fds5 2> "%smvug%\figures\fds5.version"

call SMV_Cases.bat

pause
