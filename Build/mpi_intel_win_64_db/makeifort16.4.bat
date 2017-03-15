@echo off
set arg1=%1

:: setup compiler environment
if x%arg1% == xbot goto skip1
REM call ..\..\Utilities\Scripts\setup_intel_compilers.bat
call C:\GIT\01_FDS\scripts\ifort16.4-intel64-compilervars.bat
:skip1

Title Building debug FDS (mpi) for 64 bit Windows

make SHELL="%ComSpec%" VPATH="../../Source" -f ..\makefile mpi_intel_win_64_db
if x%arg1% == xbot goto skip2
pause
:skip2
