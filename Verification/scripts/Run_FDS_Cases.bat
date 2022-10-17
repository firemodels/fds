@echo off

set curdir=%CD%
set DEBUG=
set rundebug=0
set run_all=1
set run_subset=0

set stopscript=0
call :getopts %*
cd %curdir%
if %stopscript% == 1 (
  exit /b
)


:: setup environment variables

set CURDIR="%CD%"
cd ..
set BASEDIR="%CD%"
cd ..\..
set "SVNROOT=%CD%"
cd %BASEDIR%\scripts
set SCRIPT_DIR="%CD%"
cd %BASEDIR%
set "TIME_FILE=%BASEDIR%\fds_case_times.txt"
set "waitfile=%BASEDIR%\waitfile.txt"

::*** uncomment following two lines to use OpenMP

:: set OMP_NUM_THREADS=1

:: set up environment variables for making and erasing stop files and running fds and cfast

set RUNFDS_R=call %SVNROOT%\fds\Utilities\Scripts\runfds.bat
set RUNTFDS_R=call %SVNROOT%\fds\Utilities\Scripts\runfds.bat

set RUNFDS_M=call %SVNROOT%\fds\Verification\scripts\make_stop.bat
set RUNTFDS_M=call %SVNROOT%\fds\Verification\scripts\make_stop.bat

set RUNFDS_E=call %SVNROOT%\fds\Verification\scripts\erase_stop.bat
set RUNTFDS_E=call %SVNROOT%\fds\Verification\scripts\erase_stop.bat

:: program locations

::set SH2BAT=%SVNROOT%\smv\Build\sh2bat\intel_win_64\sh2bat
set SH2BAT=sh2bat

set FDSBASE=fds_impi_intel_win%DEBUG%.exe
set "FDSEXE=%SVNROOT%\fds\Build\impi_intel_win%DEBUG%\%FDSBASE%"

:: ---------- Ensure that various programs exists

call :is_file_installed %FDSEXE%|| exit /b 1
call :is_file_installed %sh2bat%|| exit /b 1

set "FDS=%FDSEXE% "

echo.
echo Creating verification case list

if "%run_all%" == "1" (
  cd %BASEDIR%
  %SH2BAT% FDS_Cases.sh FDS_Cases.bat
)
if "%run_subset%" == "1" (
  cd %BASEDIR%
  %SH2BAT% FDS_Cases_Subset.sh FDS_Cases_Subset.bat
)

if "%rundebug%" == "1" (
  SET QFDS=%RUNFDS_M%
  SET RUNTFDS=%RUNTFDS_M%
) else (
  SET QFDS=%RUNFDS_E%
  SET RUNTFDS=%RUNTFDS_E%
)

:: create or erase stop files

if "%run_all%" == "1" (
  cd %BASEDIR%
  call FDS_Cases.bat
)
if "%run_subset%" == "1" (
  cd %BASEDIR%
  call FDS_Cases_Subset.bat
)

echo "test cases begin" >> %TIME_FILE%
date /t >> %TIME_FILE%
time /t >> %TIME_FILE%

SET QFDS=%RUNFDS_R%
SET RUNTFDS=%RUNTFDS_R%

set MPIEXEC_PORT_RANGE=
set MPICH_PORT_RANGE=

call :stop_prog hydra_bstrap_proxy.exe
call :stop_prog mpiexec.exe
call :stop_prog %FDSBASE%

if "%run_all%" == "1" (
  echo.
  echo Running FDS cases
  echo.
  cd %BASEDIR%
  call FDS_Cases.bat
)

if "%run_subset%" == "1" (
  echo.
  echo Running a Subset of FDS cases
  echo.
  cd %BASEDIR%
  call FDS_Cases_Subset.bat
)

call :wait_until_finished
goto eof

:: -------------------------------------------------------------
:stop_prog
:: -------------------------------------------------------------
set prog=%1
tasklist | grep %prog%  | wc -l > items.txt
set /p nitems=<items.txt
if %nitems% == 0 goto skip_kill
  echo stopping %prog%
  taskkill /f /im %prog% >Nul
:skip_kill
erase items.txt
exit /b

:: -------------------------------------------------------------
:wait_until_finished
:: -------------------------------------------------------------
:loop1
:: FDSBASE defined in Run_FDS_Cases
tasklist | find /i /c "%FDSBASE%" > %waitfile%
set /p numexe=<%waitfile%
echo Number of cases running - %numexe%
if %numexe% == 0 goto finished
Timeout /t 30 >nul
goto loop1
erase %waitfile%

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

:getopts
 if (%1)==() exit /b
 set valid=0
 set arg=%1
 if /I "%1" EQU "-help" (
   call :usage
   set stopscript=1
   exit /b
 )
 if /I "%1" EQU "-subset" (
   set valid=1
   set run_all=0
   set run_subset=1
 )
 if /I "%1" EQU "-debug" (
   set valid=1
   set rundebug=1
   set DEBUG=_db
 )
 shift
 if %valid% == 0 (
   echo.
   echo ***Error: the input argument %arg% is invalid
   echo.
   echo Usage:
   call :usage
   set stopscript=1
   exit /b
 )
if not (%1)==() goto getopts
exit /b

:usage
echo Run_FDS_Cases [options]
echo.
echo -help   - display this message
echo -debug  - run cases using the debug version of FDS
echo -subset - run a subset of FDS cases
exit /b

:eof
echo "FDS test cases end" >> %TIME_FILE%
date /t >> %TIME_FILE%
time /t >> %TIME_FILE%

:eof2


