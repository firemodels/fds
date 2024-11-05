@echo off
set arg1=%1

:: build hypre and/or sundials libraries if necessary

call ..\Scripts\build_thirdparty_libs %*
if %stopscript% == 1 exit /b

for %%I in (.) do set TARGET=%%~nxI

:: setup compiler environment
if x%arg1% == xbot goto endif1
call ..\Scripts\setup_intel_compilers.bat
:endif1

Title Building OpenMP FDS (Intel MPI/%INTEL_IFORT% OpenMP) for 64 bit Windows %TARGET%

make SHELL="%ComSpec%" VPATH="../../Source" -f ..\makefile %TARGET%
if x%arg1% == xbot goto endif2
pause
:endif2
