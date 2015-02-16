@echo off

:: $Date$ 
:: $Revision$
:: $Author$

set rundebug=%1
if "%rundebug%" == "1" (
  set DEBUG=_db
) else (
  set DEBUG=
)

:: setup environment variables

set BASEDIR="%CD%"
cd ..
set SVNROOT="%CD%"
cd %BASEDIR%\scripts
set SCRIPT_DIR="%CD%"
cd %SVNROOT%\..\cfast\
set CFAST=%CD%
cd %BASEDIR%
set TIME_FILE="%BASEDIR%\fds_case_times.txt"

::*** uncomment following two lines to use OpenMP

:: set OMP_NUM_THREADS=1

:: set up environment variables for making and erasing stop files and running fds and cfast

set RUNFDS_R=call %SVNROOT%\Utilities\Scripts\runfds.bat
set RUNTFDS_R=call %SVNROOT%\Utilities\Scripts\runfds.bat
set RUNCFAST_R=call %SVNROOT%\Utilities\Scripts\runcfast.bat

set RUNFDS_M=call %SVNROOT%\Verification\scripts\make_stop.bat
set RUNTFDS_M=call %SVNROOT%\Verification\scripts\make_stop.bat
set RUNCFAST_M=call %SVNROOT%\Verification\scripts\make_stop.bat

set RUNFDS_E=call %SVNROOT%\Verification\scripts\erase_stop.bat
set RUNTFDS_E=call %SVNROOT%\Verification\scripts\erase_stop.bat
set RUNCFAST_E=call %SVNROOT%\Verification\scripts\erase_stop.bat

:: program locations

set BACKGROUNDEXE=%SVNROOT%\Utilities\background\intel_win_32\background.exe
set FDSBASE=fds_mpi_win_64%DEBUG%.exe
set FDSEXE=%SVNROOT%\FDS_Compilation\mpi_intel_win_64%DEBUG%\%FDSBASE%
set CFASTEXE=%CFAST%\CFAST\intel_win_64\cfast6_win_64.exe
set WIND2FDSEXE=%SVNROOT%\Utilities\wind2fds\intel_win_64\wind2fds_win_64.exe
set SH2BAT=%SVNROOT%\Utilities\Data_Processing\sh2bat

:: ---------- Ensure that various programs exists

call :is_file_installed %BACKGROUNDEXE%|| exit /b 1
call :is_file_installed %FDSEXE%|| exit /b 1
call :is_file_installed %CFASTEXE%|| exit /b 1
call :is_file_installed %WIND2FDSEXE%|| exit /b 1
call :is_file_installed %sh2bat%|| exit /b 1

set "bg=%BACKGROUNDEXE% -u 60 -m 70 -d 5 "
set FDS=%bg%%FDSEXE%
set CFAST=%bg%%CFASTEXE%

echo.
echo Creating FDS verification case list
cd %BASEDIR%
%SH2BAT% FDS_Cases.sh FDS_Cases.bat
cd %SCRIPT_DIR%
echo Creating smokeview verification case list
%SH2BAT% SMV_Cases.sh SMV_Cases.bat
%SH2BAT% SMV_geom_Cases.sh SMV_geom_Cases.bat

echo Converting wind data
echo .
cd %SVNROOT%\Verification\WUI
%WIND2FDSEXE% -prefix sd11 -offset " 50.0  50.0 0.0" wind_data1a.csv
%WIND2FDSEXE% -prefix sd12 -offset " 50.0 150.0 0.0" wind_data1b.csv
%WIND2FDSEXE% -prefix sd21 -offset "150.0  50.0 0.0" wind_data1c.csv
%WIND2FDSEXE% -prefix sd22 -offset "150.0 150.0 0.0" wind_data1d.csv

set smvug="%SVNROOT%\Manuals\SMV_User_Guide\"
echo | %FDSEXE% 2> "%smvug%\SCRIPT_FIGURES\fds.version"

if "%rundebug%" == "1" (
  SET QFDS=%RUNFDS_M%
  SET RUNTFDS=%RUNTFDS_M%
  SET RUNCFAST=%RUNCFAST_M%
) else (
  SET QFDS=%RUNFDS_E%
  SET RUNTFDS=%RUNTFDS_E%
  SET RUNCFAST=%RUNCFAST_E%
)

:: create or erase stop files

cd %BASEDIR%
call FDS_Cases.bat
call %SCRIPT_DIR%\SMV_Cases.bat
call %SCRIPT_DIR%\SMV_geom_Cases.bat
call :wait_until_finished

echo.
echo Running FDS cases
echo.

echo "FDS test cases begin" >> %TIME_FILE%
date /t >> %TIME_FILE%
time /t >> %TIME_FILE%

SET QFDS=%RUNFDS_R%
SET RUNTFDS=%RUNTFDS_R%
SET RUNCFAST=%RUNCFAST_R%

cd %BASEDIR%
call FDS_Cases.bat
call %SCRIPT_DIR%\SMV_Cases.bat
call %SCRIPT_DIR%\SMV_geom_Cases.bat

goto eof

:: -------------------------------------------------------------
:wait_until_finished
:: -------------------------------------------------------------
:loop1
:: FDSBASE defined in Run_SMV_Cases and Run_FDS_Cases (the same in each)
tasklist | find /i /c "%FDSBASE%" > %waitfile%
set /p numexe=<%waitfile%
echo Number of cases running - %numexe%
if %numexe% == 0 goto finished
Timeout /t 30 >nul 
goto loop1

:finished
exit /b

:: -----------------------------------------
:is_file_installed
:: -----------------------------------------
  set program=%1
  %program% -help 1> %BASEDIR%\exist.txt 2>&1
  type %BASEDIR%\exist.txt | find /i /c "not recognized" > %BASEDIR%\count.txt
  set /p nothave=<%BASEDIR%\count.txt
  if %nothave% GTR 0 (
    echo "***Fatal error: %program% not present"
    echo "Verification suite aborted"
    exit /b 1
  )
  exit /b 0


:eof
echo "FDS test cases end" >> %TIME_FILE%
date /t >> %TIME_FILE%
time /t >> %TIME_FILE%

:eof2


