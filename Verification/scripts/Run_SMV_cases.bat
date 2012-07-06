@echo off
set svn_drive=c:

set SCRIPT_DIR=%CD%
set BASEDIR=%CD%\..
set SVNROOT=%BASEDIR%\..\
set TIME_FILE=%SCRIPT_DIR%\smv_case_times.txt

set RUNFDS=call %SVNROOT%\Utilities\Scripts\runfds_win32.bat
set RUNCFAST=call %SVNROOT%\Utilities\Scripts\runcfast_win32.bat

Rem VVVVVVVVVVVV set parameters VVVVVVVVVVVVVVVVVVVVVV

Rem Choose FDS version (repository or release)

set FDSEXE=%SVNROOT%\FDS_Compilation\intel_win_64\fds_win_64
Rem set FDSEXE=fds

Rem Choose CFAST version (repository or release)

set CFASTEXE=cfast6
Rem set CFASTEXE=%SVNROOT%\..\cfast\CFAST\intel_win_64\cfast6_win_64

Rem Run jobs in background (or not)

set background=background -u 85 -d 5
Rem set background=

Rem ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

set FDS=%background% %FDSEXE%
set CFAST=%background% %CFASTEXE%

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

cd %BASEDIR%
echo "smokeview test cases end" >> %TIME_FILE%
date /t >> %TIME_FILE%
time /t >> %TIME_FILE%

echo "FDS/CFAST cases completed"

pause
