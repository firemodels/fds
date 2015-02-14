@echo off

set rundebug=%1
set runonlygeom=%2

set svn_drive=c:
if "%rundebug%" == "1" (
set DEBUG=_db
) else (
set DEBUG=
)

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

set RUNFDS_R=call %SVNROOT%\Utilities\Scripts\runfds.bat
set RUNWFDS_R=call %SVNROOT%\Utilities\Scripts\runwfds.bat
set RUNTFDS_R=call %SVNROOT%\Utilities\Scripts\runfds.bat
set RUNCFAST_R=call %SVNROOT%\Utilities\Scripts\runcfast.bat

set RUNFDS_M=call %SVNROOT%\Verification\scripts\make_stop.bat
set RUNWFDS_M=call %SVNROOT%\Verification\scripts\make_stop.bat
set RUNTFDS_M=call %SVNROOT%\Verification\scripts\make_stop.bat
set RUNCFAST_M=call %SVNROOT%\Verification\scripts\make_stop.bat

set RUNFDS_E=call %SVNROOT%\Verification\scripts\erase_stop.bat
set RUNWFDS_E=call %SVNROOT%\Verification\scripts\erase_stop.bat
set RUNTFDS_E=call %SVNROOT%\Verification\scripts\erase_stop.bat
set RUNCFAST_E=call %SVNROOT%\Verification\scripts\erase_stop.bat

:: VVVVVVVVVVVV set parameters VVVVVVVVVVVVVVVVVVVVVV

set FDSBASE=fds_mpi_win_64%DEBUG%.exe
set FDSEXE=%SVNROOT%\FDS_Compilation\mpi_intel_win_64%DEBUG%\%FDSBASE%
set CFASTEXE=%CFAST%\CFAST\intel_win_64\cfast6_win_64.exe
set WIND2FDSEXE=%SVNROOT%\Utilities\wind2fds\intel_win_64\wind2fds_win_64.exe

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
set CFAST=%bg%%CFASTEXE%

set SH2BAT=%SVNROOT%\Utilities\Data_Processing\sh2bat
call :is_file_installed %sh2bat%|| exit /b 1

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
%SH2BAT% SMV_geom_Cases.sh SMV_geom_Cases.bat

cd %BASEDIR%

echo "smokeview test cases begin" > %TIME_FILE%
date /t >> %TIME_FILE%
time /t >> %TIME_FILE%

:: create a text file containing the FDS version used to run these tests.
:: This file is included in the smokeview user's guide

set smvug="%SVNROOT%\Manuals\SMV_User_Guide\"
echo | %FDSEXE% 2> "%smvug%\SCRIPT_FIGURES\fds.version"

if "%rundebug%" == "1" (
  SET QFDS=%RUNFDS_M%
  SET RUNWFDS=%RUNWFDS_M%
  SET RUNTFDS=%RUNTFDS_M%
  SET RUNCFAST=%RUNCFAST_M%
) else (
  SET QFDS=%RUNFDS_E%
  SET RUNWFDS=%RUNWFDS_E%
  SET RUNTFDS=%RUNTFDS_E%
  SET RUNCFAST=%RUNCFAST_E%
)

:: create or erase stop files

if "%runonlygeom%" == "1" (
  call %SCRIPT_DIR%\SMV_geom_Cases.bat
) else (
  call %SCRIPT_DIR%\SMV_Cases.bat
  call %SCRIPT_DIR%\SMV_geom_Cases.bat
)

:: run cases

SET QFDS=%RUNFDS_R%
SET RUNWFDS=%RUNWFDS_R%
SET RUNTFDS=%RUNTFDS_R%
SET RUNCFAST=%RUNCFAST_R%
if "%runonlygeom%" == "1" (
  call %SCRIPT_DIR%\SMV_geom_Cases.bat
) else (
  call %SCRIPT_DIR%\SMV_Cases.bat
  call %SCRIPT_DIR%\SMV_geom_Cases.bat
)

cd %BASEDIR%
echo "smokeview test cases end" >> %TIME_FILE%
date /t >> %TIME_FILE%
time /t >> %TIME_FILE%

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
