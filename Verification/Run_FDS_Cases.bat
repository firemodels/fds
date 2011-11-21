@echo off
set svn_drive=d:

set BASEDIR="%CD%"
set SVNROOT=%BASEDIR%\..\
set TIME_FILE="%BASEDIR%\fds_case_times.txt"

Rem Choose one of the following four FDS "definitions" by commenting all lines but one.

set FDSEXE=%SVNROOT%\FDS_Compilation\intel_win_32\fds_win_32
Rem set FDSEXE=%SVNROOT%\FDS_Compilation\intel_win_64\fds_win_64


Rem Choose one of the following run options by commenting the line you don't want to use.
Rem Use the first "FDS" definition to run one case at a time in the forground.
Rem Use the second "FDS" definition to execute more than one case at a time in the background

set FDS=%FDSEXE%
Rem set FDS=background -u 75 -d 10 %FDSEXE%


set RUNFDS=call %SVNROOT%\Utilities\Scripts\runfds_win32.bat

echo You are about to run the Verification Test Suite.
echo Press any key to begin.
pause > Nul

Rem -------------------------
Rem run FDS veriication cases
Rem -------------------------
echo creating FDS case list from FDS_Cases.csh
..\Utilities\Data_Processing\sh2bat FDS_Cases.sh FDS_Cases.bat

echo running FDS cases

echo "FDS test cases begin" > %TIME_FILE%
date /t >> %TIME_FILE%
time /t >> %TIME_FILE%

call FDS_Cases.bat

echo "FDS test cases end" > %TIME_FILE%
date /t >> %TIME_FILE%
time /t >> %TIME_FILE%

erase FDS_Cases.bat

pause
