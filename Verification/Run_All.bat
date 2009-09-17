@echo off
set svn_drive=d:
set SVNROOT=d:\FDS-SMV
set FDS=%SVNROOT%\FDS_Compilation\intel_win_32\fds5_win_32
set RUNFDS=call %SVNROOT%\Utilities\Scripts\runfds_win32.bat
set BASEDIR=%CD%
%svn_drive%

echo You are about to run the Verification Test Suite.
echo Press any key begin.
pause > Nul

call ./FDS_Cases.bat

echo FDS cases submitted
