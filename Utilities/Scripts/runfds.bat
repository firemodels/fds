@echo off

:: set number of openmp threads

set OMP_NUM_THREADS=1

set scriptdir=%~p0
set CURDIR=%CD%
cd %scriptdir%\..\..\..
set ROOT=%CD%
cd %CURDIR%
set njob_file=njob_file
set MAX_CPU_LOAD=10

call %ROOT%\fds\Utilities\Scripts\getopts.bat %*

set fulldir=%BASEDIR%/%dir%
For %%A in ("%FDS%") do (
   set PROG=%%~nxA
)

set in=%infile%.fds
set stopfile=%infile%.stop
set runfds_start=%ROOT%\fds\Utilities\Scripts\runfds_start

cd %fulldir%


if exist %stopfile% (
   erase %stopfile%
)
if "%rundebug%" == "1" (
   echo 2 > %stopfile%
)

call :pause_fds

echo %in% started
echo.
start /min "%infile%" call %runfds_start% %FDS% %in% %infile%.err

goto eof

:: -------------------------------------------------------------
:pause_fds
:: -------------------------------------------------------------

call :get_cpu_load
if %loadavg% LSS %MAX_CPU_LOAD% goto exit_pause_fds

:loop2
echo starting %in% when CPU load ^< %MAX_CPU_LOAD%
:loop1

call :get_cpu_load
if %loadavg% LSS %MAX_CPU_LOAD% goto exit_pause_fds

Timeout /t 10 >nul
goto loop1

:exit_pause_fds
exit /b

:: -------------------------------------------------------------
:get_cpu_load
:: -------------------------------------------------------------

:: FDSBASE defined in Run_FDS_Cases
tasklist | find /i /c "%FDSBASE%" > %waitfile%
set /p loadavg=<%waitfile%
erase %waitfile%
echo Number of cases running - %loadavg%
Timeout /t 6 >nul

exit /b

:eof

