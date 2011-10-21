@echo off

Rem windows batch file to build smokeview from the command line

IF "%SETUP_IFORT_COMPILER12%"=="1" GOTO envexist

set SETUP_IFORT_COMPILER12=1

echo Setting up compiler environment
call "%IFORT_COMPILER12%\bin\ifortvars" ia32_intel64
:envexist
erase *.obj
erase *.mod
make VPATH="../../FDS_Source" -f ..\makefile intel_win_64
pause