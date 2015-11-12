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

set CURDIR="%CD%"
cd ..
set BASEDIR="%CD%"
cd ..
set SVNROOT="%CD%"
cd %BASEDIR%\scripts
set SCRIPT_DIR="%CD%"
cd %BASEDIR%
set TIME_FILE="%BASEDIR%\fds_case_times.txt"

::*** uncomment following two lines to use OpenMP

:: set OMP_NUM_THREADS=1

:: set up environment variables for making and erasing stop files and running fds and cfast

set RUNFDS_R=call %SVNROOT%\Utilities\Scripts\runfds.bat
set RUNTFDS_R=call %SVNROOT%\Utilities\Scripts\runfds.bat

set RUNFDS_M=call %SVNROOT%\Verification\scripts\make_stop.bat
set RUNTFDS_M=call %SVNROOT%\Verification\scripts\make_stop.bat

set RUNFDS_E=call %SVNROOT%\Verification\scripts\erase_stop.bat
set RUNTFDS_E=call %SVNROOT%\Verification\scripts\erase_stop.bat

:: program locations

set BACKGROUNDEXE=%SVNROOT%\Utilities\background\intel_win_64\background.exe
set FDSBASE=fds_mpi_win_64%DEBUG%.exe
set FDSEXE=%SVNROOT%\FDS_Compilation\mpi_intel_win_64%DEBUG%\%FDSBASE%
set SH2BAT=%SVNROOT%\Utilities\Data_Processing\sh2bat

:: ---------- Ensure that various programs exists

call :is_file_installed %BACKGROUNDEXE%|| exit /b 1
call :is_file_installed %FDSEXE%|| exit /b 1
call :is_file_installed %sh2bat%|| exit /b 1

set "bg=%BACKGROUNDEXE% -u 60 -m 70 -d 5 "
set FDS=%bg%%FDSEXE%

echo.
echo Creating FDS verification case list
cd %BASEDIR%
%SH2BAT% FDS_Cases.sh FDS_Cases.bat

if "%rundebug%" == "1" (
  SET QFDS=%RUNFDS_M%
  SET RUNTFDS=%RUNTFDS_M%
) else (
  SET QFDS=%RUNFDS_E%
  SET RUNTFDS=%RUNTFDS_E%
)

:: create or erase stop files

cd %BASEDIR%
call FDS_Cases.bat

echo.
echo Running FDS cases
echo.

echo "FDS test cases begin" >> %TIME_FILE%
date /t >> %TIME_FILE%
time /t >> %TIME_FILE%

SET QFDS=%RUNFDS_R%
SET RUNTFDS=%RUNTFDS_R%

cd %BASEDIR%
call FDS_Cases.bat
call :wait_until_finished

goto eof

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


