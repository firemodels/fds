@echo off
set svn_drive=c:

set SCRIPT_DIR=%CD%
set BASEDIR=%CD%\..
set SVNROOT=%BASEDIR%\..\
set TIME_FILE=%SCRIPT_DIR%\smv_case_times.txt

Rem Choose one of the following four FDS "definitions" by commenting all lines but one.

Rem set FDSEXE=%SVNROOT%\FDS_Compilation\intel_win_32\fds_win_32
Rem set FDSEXE=fds_win_32
set FDSEXE=%SVNROOT%\FDS_Compilation\intel_win_64\fds_win_64
Rem set FDSEXE=fds5_win_64

Rem Choose one of the following run options by commenting the line you don't want to use.
Rem Use the first "FDS" definition to run one case at a time in the forground.
Rem Use the second "FDS" definition to execute more than one case at a time in the background

Rem set FDS=%FDSEXE%
set FDS=background -u 85 -d 5 %FDSEXE%

set RUNFDS=call %SVNROOT%\Utilities\Scripts\runfds_win32.bat

echo You are about to run the Smokeview Verification Test Suite.
echo Press any key to proceed, CTRL c to abort
pause > Nul

echo creating FDS case list from SMV_Cases.sh
..\..\Utilities\Data_Processing\sh2bat SMV_Cases.sh SMV_Cases.bat

cd %BASEDIR%

echo "smokeview test cases begin" > %TIME_FILE%
date /t >> %TIME_FILE%
time /t >> %TIME_FILE%

Rem create a text file containing the FDS5 version used to run these tests.
Rem This file is included in the smokeview user's guide

set smvug="%CD%\..\Manuals\SMV_User_Guide\"
echo | %FDSEXE% 2> "%smvug%\figures\fds5.version"

call %SCRIPT_DIR%\SMV_Cases.bat

erase %SCRIPT_DIR%\SMV_Cases.bat

cd %SCRIPT_DIR%
call run_wui_tree_test

cd %BASEDIR%
echo "smokeview test cases end" >> %TIME_FILE%
date /t >> %TIME_FILE%
time /t >> %TIME_FILE%

echo "FDS cases completed"

pause
