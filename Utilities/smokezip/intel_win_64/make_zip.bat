@echo off

:: windows batch file to build smokezip from the command line

IF "%SETUP_IFORT_COMPILER%"=="1" GOTO envexist

set SETUP_IFORT_COMPILER=1

echo Setting up compiler environment
call "%IFORT_COMPILER14%\bin\compilervars" intel64
:envexist

Title Building smokezip for 64 bit Windows
erase *.obj *.mod
make -f ..\Makefile intel_win_64
pause+

