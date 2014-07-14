@echo off

:: windows batch file to build wind2fds from the command line

IF "%VS_VERSION%"=="" SET VS_VERSION=vs2012
IF "%SETUP_IFORT_COMPILER32%"=="1" GOTO envexist

set SETUP_IFORT_COMPILER32=1

echo Setting up compiler environment
call "%IFORT_COMPILER14%\bin\compilervars" ia32
:envexist
Title Building wind2fds for 32 bit Windows
erase *.obj
erase *.mod
make -f ..\Makefile intel_win_32

