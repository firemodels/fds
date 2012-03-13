@echo off

Rem windows batch file to build smokeview from the command line

IF "%SETUP_IFORT_COMPILER_IA32%"=="1" GOTO envexist

set SETUP_IFORT_COMPILER_IA32=1

echo Setting up compiler environment
call "%IFORT_COMPILER12%\bin\ifortvars" intel64 vs2008
call "%IFORT_COMPILER12%\bin\iclvars" intel64 vs2008
call "%IFORT_COMPILER12%\bin\ifortvars" ia32 vs2008
call "%IFORT_COMPILER12%\bin\iclvars" ia32 vs2008
:envexist
make -f ..\Makefile intel_win_32
pause