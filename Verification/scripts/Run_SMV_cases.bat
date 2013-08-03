@echo off
set svn_drive=c:

set SCRIPT_DIR=%CD%
cd %CD%\..
set BASEDIR=%CD%

cd %BASEDIR%\..\
set SVNROOT=%CD%

cd %SVNROOT%\..\cfast\
set CFAST=%CD%

set TIME_FILE=%SCRIPT_DIR%\smv_case_times.txt

set RUNFDS=call %SVNROOT%\Utilities\Scripts\runfds_win32.bat
set RUNCFAST=call %SVNROOT%\Utilities\Scripts\runcfast_win32.bat

Rem VVVVVVVVVVVV set parameters VVVVVVVVVVVVVVVVVVVVVV

Rem Choose FDS version (repository or release)

set FDSEXE=%SVNROOT%\FDS_Compilation\intel_win_64\fds_win_64.exe
Rem set FDSEXE=%SVNROOT%\FDS_Compilation\intel_win_64_db\fds_win_64_db.exe
Rem set FDSEXE=fds

set BACKGROUNDEXE=%SVNROOT%\Utilities\background\intel_win_32\background.exe

Rem Choose CFAST version (repository or release)

Rem set CFASTEXE=cfast6
set CFASTEXE=%CFAST%\CFAST\intel_win_64\cfast6_win_64.exe


Rem Define wind2fds

set WIND2FDSEXE=%SVNROOT%\Utilities\wind2fds\intel_win_64\wind2fds_win_64.exe

Rem Run jobs in background (or not)

set "bg=%BACKGROUNDEXE% -u 85 -d 5 "
Rem set bg=

Rem ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Rem ---------- Ensure that wind2fds_win_64 exists

IF EXIST %WIND2FDSEXE% GOTO endif_wind2fdsexist
echo ***Fatal error.  The file %WIND2FDSEXE% does not exist. 
echo Aborting now...
pause>NUL
goto:eof

:endif_wind2fdsexist

Rem ---------- Ensure that fds exists

IF EXIST %FDSEXE% GOTO endif_fdsexist
echo ***Fatal error.  The file %FDSEXE% does not exist. 
echo Aborting now...
pause>NUL
goto:eof

:endif_fdsexist

Rem ---------- Ensure that CFAST exists

IF EXIST %CFASTEXE% GOTO endif_cfastexist
echo ***Fatal error.  The file %CFASTEXE% does not exist. 
echo Aborting now...
pause>NUL
goto:eof

:endif_cfastexist

set FDS=%bg%%FDSEXE%
set CFAST=%bg%%CFASTEXE%
set SH2BAT=%SVNROOT%\Utilities\Data_Processing\sh2bat

echo You are about to run the Smokeview Verification Test Suite.
echo.
echo FDS=%FDS%
echo CFAST=%CFAST%
echo.
echo Press any key to proceed, CTRL c to abort
pause > Nul

echo Converting wind data
echo .
cd %SVNROOT%\Verification\WUI
%WIND2FDSEXE% -prefix sd11 -offset " 50.0  50.0 0.0" wind_data1a.csv
%WIND2FDSEXE% -prefix sd12 -offset " 50.0 150.0 0.0" wind_data1b.csv
%WIND2FDSEXE% -prefix sd21 -offset "150.0  50.0 0.0" wind_data1c.csv
%WIND2FDSEXE% -prefix sd22 -offset "150.0 150.0 0.0" wind_data1d.csv

cd %SCRIPT_DIR%
echo creating FDS case list from SMV_Cases.sh
%SH2BAT% SMV_Cases.sh SMV_Cases.bat

cd %BASEDIR%

echo "smokeview test cases begin" > %TIME_FILE%
date /t >> %TIME_FILE%
time /t >> %TIME_FILE%

Rem create a text file containing the FDS5 version used to run these tests.
Rem This file is included in the smokeview user's guide

set smvug="%SVNROOT%\Manuals\SMV_User_Guide\"
echo | %FDSEXE% 2> "%smvug%\SCRIPT_FIGURES\fds.version"

call %SCRIPT_DIR%\SMV_Cases.bat

erase %SCRIPT_DIR%\SMV_Cases.bat

cd %BASEDIR%
echo "smokeview test cases end" >> %TIME_FILE%
date /t >> %TIME_FILE%
time /t >> %TIME_FILE%

echo "FDS/CFAST cases completed"

:eof
pause
