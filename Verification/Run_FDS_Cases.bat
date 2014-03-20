@echo off

set OPENMP=

:: to use openmp (assuming you have compiled the openmp exe's), uncomment the follwoing two "set" lines
:: and set the number of threads to what you want (or leave commented to use number of cores available)
:: set OPENMP=openmp_
:: set OMP_NUM_THREADS=2

set BASEDIR="%CD%"
set SVNROOT=%BASEDIR%\..\
set TIME_FILE="%BASEDIR%\fds_case_times.txt"

:: VVVVVVVVVVVVVVVVVVVVVV define FDS exe location VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

:: DO NOT erase or comment the following two line when committing to repository

set FDSBASE=fds_%OPENMP%win_64.exe
set FDSEXE=%SVNROOT%\FDS_Compilation\%OPENMP%intel_win_64\fds_%OPENMP%win_64.exe

:: if you wish to use a custom FDS, comment the above two "set" lines and uncomment
:: and customize the below two "set" lines

:: set FDSBASE=fds_win_64.exe
:: set FDSEXE=C:\PROJECTS\FDS_DEVEL\6.0a\fds_win_64

:: ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

set BACKGROUNDEXE=%SVNROOT%\Utilities\background\intel_win_32\background.exe

Rem Choose one of the following run options by commenting the line you don't want to use.
Rem Use the first "FDS" definition to run one case at a time in the forground.
Rem Use the second "FDS" definition to execute more than one case at a time in the background

Rem set FDS=%FDSEXE%
set FDS=%BACKGROUNDEXE% -u 75 -d 10 %FDSEXE%

set RUNFDS=call %SVNROOT%\Utilities\Scripts\runfds_win32.bat
set RUNFDSMPI=call %SVNROOT%\Utilities\Scripts\runfdsmpi_win32.bat

echo You are about to run the Verification Test Suite.
echo Press any key to begin.
pause > Nul

echo.
echo Creating FDS case list from FDS_Cases.sh
..\Utilities\Data_processing\sh2bat FDS_Cases.sh FDS_Cases.bat
echo Creating FDS_MPI case list from FDS_MPI_Cases.sh
..\Utilities\Data_processing\sh2bat FDS_MPI_Cases.sh FDS_MPI_Cases.bat

echo.
echo Running FDS cases
echo.

echo "FDS test cases begin" >> %TIME_FILE%
date /t >> %TIME_FILE%
time /t >> %TIME_FILE%

call FDS_Cases.bat

call FDS_MPI_Cases.bat

:: loop until all FDS cases have finished

:loop1
tasklist | find /i /c "%FDSBASE%" > temp.out
set /p numexe=<temp.out
echo Number of cases running - %numexe%
if %numexe% == 0 goto finished
Timeout /t 30 >nul 
goto loop1

:finished
echo "FDS cases completed"
goto eof

::  for possible future edits, put any batch subroutines added between here and eof statement below

:eof
echo "FDS test cases end" >> %TIME_FILE%
date /t >> %TIME_FILE%
time /t >> %TIME_FILE%

erase FDS_Cases.bat
erase FDS_MPI_Cases.bat

pause
