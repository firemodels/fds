@echo off
set intelbin="%IFORT_COMPILER14%\bin"

IF "%SETUP_IFORT_COMPILER_64%"=="1" GOTO envexist

set SETUP_IFORT_COMPILER_64=1

echo Setting up compiler environment
call %intelbin%\ifortvars intel64

:envexist
Title Building FDS for 64 bit Windows
make VPATH="../../FDS_Source" -f ..\makefile intel_win_64
pause

