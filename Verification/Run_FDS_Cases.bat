@echo off
set svn_drive=d:

set BASEDIR="%CD%"
set SVNROOT=%BASEDIR%\..\
set TIME_FILE="%BASEDIR%\fds_case_times.txt"

Rem Choose one of the following four FDS "definitions" by commenting all lines but one.

set FDSEXE=\projects\fds_devel\6.0a\fds_win_64
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

echo "FDS test cases end" >> %TIME_FILE%
date /t >> %TIME_FILE%
time /t >> %TIME_FILE%

erase FDS_Cases.bat
erase FDS_MPI_Cases.bat

pause
