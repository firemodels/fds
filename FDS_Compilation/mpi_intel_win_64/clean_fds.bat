@echo off
set arg1=%1

:: setup compiler environment
if x%arg1% == xbot goto skip1
call ..\..\Utilities\Scripts\setup_intel_compilers.bat
:skip1

Title Cleaning FDS (mpi) for 64 bit Windows

make SHELL="%ComSpec%" VPATH="../../FDS_Source" -f ..\makefile clean
if x%arg1% == xbot goto skip2
pause
:skip2
