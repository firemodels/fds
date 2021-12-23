@echo off
set arg1=%1
set arg2=%2

if x%arg1% == x goto skip2
if x%arg1% == xbot goto skip1
  set arg2=%arg1%
  set arg1=
:skip2

call ..\Scripts\set_intel_compiler.bat %arg2%

:: setup compiler environment
if x%arg1% == xbot goto skip1
call ..\Scripts\setup_intel_compilers.bat
:skip1

Title Building DV FDS (Intel MPI/%INTEL_IFORT%) for 64 bit Windows

make SHELL="%ComSpec%" VPATH="../../Source" -f ..\makefile impi_intel_win_64_dv
if x%arg1% == xbot goto skip2
pause
:skip2
