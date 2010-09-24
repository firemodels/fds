@echo off

Rem windows batch file to build background from the command line

IF "%SETUP_IFORT_COMPILER11%"=="1" GOTO envexist

set SETUP_IFORT_COMPILER11=1

echo Setting up compiler environment
call "%IFORT_COMPILER11%\bin\ifortvars" ia32
call "%IFORT_COMPILER11%\bin\iclvars" ia32

:envexist
erase *.obj
erase *.mod
make -f ..\Makefile intel_win_32
pause
