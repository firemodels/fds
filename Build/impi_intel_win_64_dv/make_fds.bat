@echo off
set arg1=%1
set arg2=%2

call ..\Scripts\set_intel_compiler.bat %arg1% %arg2%

:: setup compiler environment
if x%arg1% == xbot goto endif1
call ..\Scripts\setup_intel_compilers.bat
:endif1

Title Building DV FDS (Intel MPI/%INTEL_IFORT%) for 64 bit Windows

make SHELL="%ComSpec%" VPATH="../../Source" -f ..\makefile impi_intel_win_64_dv
if x%arg1% == xbot goto endif2
pause
:endif2
