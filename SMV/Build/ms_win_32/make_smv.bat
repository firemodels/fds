@echo off

Rem windows batch file to build smokeview from the command line

IF "%SETUP_IFORT_COMPILER_IA32%"=="1" GOTO envexist

set SETUP_IFORT_COMPILER_IA32=1

echo Setting up compiler environment
call "%IFORT_COMPILER12%\bin\ifortvars" ia32
if exist "%VS_COMPILER%\..\vcvars32.bat" call "%VS_COMPILER%\..\vcvars32.bat"
:envexist
make -j4 -f ..\Makefile ms_win_32
pause