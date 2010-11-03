@echo off

Rem windows batch file to build smokeview from the command line

IF "%SETUP_IFORT_COMPILER11%"=="1" GOTO envexist

set SETUP_IFORT_COMPILER11=1

echo Setting up compiler environment
call "%IFORT_COMPILER11%\bin\ifortvars" intel64
call "%IFORT_COMPILER11%\bin\iclvars" intel64
:envexist
erase *.obj
erase *.mod
make -j4 -f ..\Makefile intel_win_test_64
pause