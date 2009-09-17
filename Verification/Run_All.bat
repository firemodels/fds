@echo off
set svn_drive=d:
set SVNROOT=d:\FDS-SMV
set FDS=%SVNROOT%\FDS_Compilation\intel_win_32\fds5_win_32
set RUNFDS=call %SVNROOT%\Utilities\Scripts\runfds_win32.bat
set BASEDIR=%CD%
%svn_drive%

call ./FDS_Cases.bat

echo FDS cases submitted
