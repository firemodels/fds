@echo off
set svn_drive=d:

set BASEDIR=%CD%
set SVNROOT=%BASEDIR%\..\

set FDS=%SVNROOT%\FDS\Build\intel_win_32\fds5_win_32

set RUNFDS=call %SVNROOT%\Utilities\Scripts\runfds_win32.bat

echo You are about to run the Verification Test Suite for Evacuation.
pause Press any key to continue or CTRL-C to abort > Nul

echo Evacuation cases submitted

call Evacuation_Cases.bat

cd %BASEDIR%

Rem create a text file containing the FDS5 version used to run these tests.
Rem This file is included in the smokeview user's guide

REM set smvug="%CD%\..\Manuals\SMV_5_User_Guide\"
REM echo | fds5 2> "%smvug%\figures\fds5.version"

cd %BASEDIR%\Evacuation
call Make_SMV_Evac_pictures.bat

echo Evacuation Verification cases done
pause Press any key...
