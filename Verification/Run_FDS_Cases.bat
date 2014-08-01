@echo off

set BASEDIR="%CD%"
cd ..
set SVNROOT="%CD%"
cd %BASEDIR%
set TIME_FILE="%BASEDIR%\fds_case_times.txt"

::*** uncomment following two lines to use OpenMP

:: set OMP_NUM_THREADS=1

:: default FDS location

set FDSEXE=%SVNROOT%\FDS_Compilation\%intel_win_64\fds_%win_64.exe

::*** uncomment following line to define a custom FDS location - (Do not commit un-commented)

:: set FDSEXE=C:\PROJECTS\FDS_DEVEL\6.0a\fds_win_64

if not exist %FDSEXE%  (
  echo "***error: The program, %FDSEXE% , was not found.  Verification test runs aborted."
  goto eof2
)

call :getfilename %FDSEXE% 
set FDSBASE=%file%

set BACKGROUNDEXE=%SVNROOT%\Utilities\background\intel_win_32\background.exe
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

:getfilename
set file=%~nx1
exit /b

:eof
echo "FDS test cases end" >> %TIME_FILE%
date /t >> %TIME_FILE%
time /t >> %TIME_FILE%

erase FDS_Cases.bat
erase FDS_MPI_Cases.bat

:eof2

pause
