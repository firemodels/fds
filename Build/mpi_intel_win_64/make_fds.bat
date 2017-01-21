@echo off
set arg1=%1
set md5hash=..\..\Utilities\Scripts\md5hash

:: setup compiler environment
if x%arg1% == xbot goto skip1
call ..\..\Utilities\Scripts\setup_intel_compilers.bat
:skip1

Title Building FDS (mpi) for 64 bit Windows

make SHELL="%ComSpec%" VPATH="../../Source" -f ..\makefile mpi_intel_win_64
%md5hash% fds_mpi_win_64.exe
if x%arg1% == xbot goto skip2
pause
:skip2
