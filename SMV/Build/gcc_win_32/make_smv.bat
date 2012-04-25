@echo off

Rem windows batch file to build smokeview from the command line

IF "%SETUP_IFORT_COMPILER_IA32%"=="1" GOTO envexist

set SETUP_IFORT_COMPILER_IA32=1

echo Setting up compiler environment
set PATH=c:\mingw\bin:%PATH%

:envexist
make -f ..\Makefile gcc_win_32
pause