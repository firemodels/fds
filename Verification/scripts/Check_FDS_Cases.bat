@echo off

set option=%1

set full=1
set subset=0


if NOT "x%option%" == "x-subset" goto skip_subset
set full=0
set subset=1
:skip_subset

cd ..
set BASEDIR="%CD%"
cd ..\..
set SVNROOT="%CD%"
cd %BASEDIR%

set QFDS=call %SVNROOT%\fds\Verification\scripts\checkfds.bat

if "%full%" == "0" goto skip_full
  call FDS_Cases.bat
:skip_full

if "%subset%" == "0" goto skip_subset
  call FDS_Cases_Subset.bat
:skip_subset
