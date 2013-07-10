@echo off

Rem windows batch file to build smokezip from the command line

IF "%SETUP_IFORT_COMPILER64%"=="1" GOTO envexist

set SETUP_IFORT_COMPILER64=1

echo Setting up compiler environment
call "%IFORT_COMPILER13%\bin\compilervars" intel64
:envexist
erase *.obj
erase *.mod
make -f ..\Makefile intel_win_64
pause
