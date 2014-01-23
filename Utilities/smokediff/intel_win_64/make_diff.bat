@echo off

Rem windows batch file to build smokediff from the command line

IF "%SETUP_IFORT_COMPILER_INTEL64%"=="1" GOTO envexist

set SETUP_IFORT_COMPILER_INTEL64=1

echo Setting up compiler environment
call "%IFORT_COMPILER14%\bin\compilervars" intel64

:envexist
erase *.obj *.mod
make -f ..\Makefile intel_win_64
pause

