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

set FDSEXE=%SVNROOT%\FDS_Compilation\intel_win_64\fds_win_64
Rem set FDSEXE=%SVNROOT%\FDS_Compilation\intel_win_64_db\fds_win_64_db
Rem set FDSEXE=fds

Rem Choose CFAST version (repository or release)

Rem set CFASTEXE=cfast6
set CFASTEXE=%CFAST%\CFAST\intel_win_64\cfast6_win_64

Rem Run jobs in background (or not)

set "bg=background -u 85 -d 5 "
Rem set background=

Rem ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Rem Check to make fds exists

IF EXIST %FDSEXE% GOTO endif_fdsexist
echo ***Fatal error.  The file %FDSEXE% does not exist. 
echo Aborting now...
pause>NUL
goto:eof

:endif_fdsexist

Rem Check to make fds exists

IF EXIST %CFASTEXE% GOTO endif_cfastexist
echo ***Fatal error.  The file %CFASTEXE% does not exist. 
echo Aborting now...
pause>NUL
goto:eof

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
