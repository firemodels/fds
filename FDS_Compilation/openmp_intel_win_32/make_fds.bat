@echo off
set intelbin="%IFORT_COMPILER14%\bin"

IF "%SETUP_IFORT_COMPILER_32%"=="1" GOTO envexist

set SETUP_IFORT_COMPILER_32=1

echo Setting up compiler environment
call %intelbin%\ifortvars ia32

:envexist
Title Building FDS (openmp) for 32 bit Windows
make VPATH="../../FDS_Source" -f ..\makefile openmp_intel_win_32
pause
