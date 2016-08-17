@echo off

Rem windows batch file to build smokeview from the command line

IF "%SETUP_IFORT_COMPILER11%"=="1" GOTO envexist

set SETUP_IFORT_COMPILER11=1

echo Setting up compiler environment
call "%IFORT_COMPILER12%\bin\ifortvars" intel64
call "%IFORT_COMPILER12%\bin\iclvars" intel64
:envexist

erase *.obj
make -j4 -f ..\Makefile intel_win_64
pause