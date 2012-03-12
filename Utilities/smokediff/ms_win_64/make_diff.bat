@echo off

Rem windows batch file to build background from the command line

IF "%SETUP_IFORT_COMPILER12%"=="1" GOTO envexist

set SETUP_IFORT_COMPILER12=1

echo Setting up compiler environment
call "%IFORT_COMPILER12%\bin\ifortvars" intel64

:envexist
if exist "%VS_COMPILER%\vcvars32x86_amd64.bat" call "%VS_COMPILER%\vcvars32x86_amd64"
erase *.obj
erase *.mod
make -f ..\Makefile ms_win_64
pause
