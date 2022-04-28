@echo off
set arg1=%1

:: setup compiler environment
if x%arg1% == xbot goto endif1
call ..\Scripts\setup_intel_compilers.bat
:endif1

Title Building FDS (Intel MPI/%INTEL_IFORT% OpenMP) for 64 bit Windows

make SHELL="%ComSpec%" VPATH="../../Source" -f ..\makefile impi_intel_win_openmp_db
if x%arg1% == xbot goto endif2
pause
:endif2
