@echo off
set SVNROOT="C:\Projects\NIST 2821-000\SVN\FDS"
set FDS="C:\Projects\NIST 2821-000\5.0source\c2\fds532.exe"
set RUNFDS=call %SVNROOT%\Utilities\Scripts\runfds_win32.bat
set BASEDIR=%CD%

echo You are about to run the Verification Test Suite.
echo Press any key begin.
pause > Nul

call FDS_Cases.bat

echo FDS cases submitted
