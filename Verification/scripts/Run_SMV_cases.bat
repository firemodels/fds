@echo off

set size=%1

set svn_drive=c:

set SCRIPT_DIR=%CD%
cd %CD%\..
set BASEDIR=%CD%

cd %BASEDIR%\..\
set SVNROOT=%CD%

cd %SVNROOT%\..\cfast\
set CFAST=%CD%

cd %SVNROOT%\..\FIRE-LOCAL
set FIRELOCAL=%CD%

set TIME_FILE=%SCRIPT_DIR%\smv_case_times.txt

set RUNFDS=call %SVNROOT%\Utilities\Scripts\runfds_win32.bat
set RUNWFDS=call %SVNROOT%\Utilities\Scripts\runwfds_win32.bat
set RUNTFDS=call %SVNROOT%\Utilities\Scripts\runwfds_win32.bat
set RUNCFAST=call %SVNROOT%\Utilities\Scripts\runcfast_win32.bat

:: VVVVVVVVVVVV set parameters VVVVVVVVVVVVVVVVVVVVVV

:: Choose FDS version (size is "", 32 or 64)

if "%size%" == "" (
  set FDSBASE=fds.exe
  set FDSEXE=%FDSBASE%
  set CFASTEXE=cfast6
  set WIND2FDSEXE=wind2fds
) else (
  set FDSBASE=fds_win_%size%.exe
  set FDSEXE=%SVNROOT%\FDS_Compilation\intel_win_%size%\fds_win_%size%.exe
  set CFASTEXE=%CFAST%\CFAST\intel_win_%size%\cfast6_win_%size%.exe
  set WIND2FDSEXE=%SVNROOT%\Utilities\wind2fds\intel_win_%size%\wind2fds_win_%size%.exe
)

set WFDSEXE=%FIRELOCAL%\bin\wfds6_9977_win_64.exe

set BACKGROUNDEXE=%SVNROOT%\Utilities\background\intel_win_32\background.exe

:: Run jobs in background (or not)

set "bg=%BACKGROUNDEXE% -u 85 -d 5 "
:: set bg=

:: ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:: ---------- Ensure that cfast, fds and wind2fds exists

call :is_file_installed %CFASTEXE%|| exit /b 1
call :is_file_installed %FDSEXE%|| exit /b 1
call :is_file_installed %WIND2FDSEXE%|| exit /b 1

set FDS=%bg%%FDSEXE%
set WFDS=%bg%%WFDSEXE%
set CFAST=%bg%%CFASTEXE%
set SH2BAT=%SVNROOT%\Utilities\Data_Processing\sh2bat

echo.
echo FDS=%FDS%
echo WFDS=%WFDS%
echo CFAST=%CFAST%
echo.

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

:: create a text file containing the FDS version used to run these tests.
:: This file is included in the smokeview user's guide

set smvug="%SVNROOT%\Manuals\SMV_User_Guide\"
echo | %FDSEXE% 2> "%smvug%\SCRIPT_FIGURES\fds.version"

call %SCRIPT_DIR%\SMV_Cases.bat

:: erase %SCRIPT_DIR%\SMV_Cases.bat

cd %BASEDIR%
echo "smokeview test cases end" >> %TIME_FILE%
date /t >> %TIME_FILE%
time /t >> %TIME_FILE%

:loop1
tasklist | find /i /c "%FDSBASE%" > temp.out
set /p numexe=<temp.out
echo Number of cases running - %numexe%
if %numexe% == 0 goto finished
Timeout /t 30 >nul 
goto loop1

:finished
echo "FDS/CFAST cases completed"
goto eof

:: -----------------------------------------
:is_file_installed
:: -----------------------------------------
  
  set program=%1
  %program% -help 1> %SCRIPT_DIR%\exist.txt 2>&1
  type %SCRIPT_DIR%\exist.txt | find /i /c "not recognized" > %SCRIPT_DIR%\count.txt
  set /p nothave=<%SCRIPT_DIR%\count.txt
  if %nothave% GTR 0 (
    echo "***Fatal error: %program% not present"
    echo "Verification suite aborted"
    exit /b 1
  )
  exit /b 0

:eof
