@echo off

Rem windows batch file to build smokeview from the command line

IF "%SETUP_IFORT_COMPILER_64%"=="1" GOTO envexist

set SETUP_IFORT_COMPILER_64=1

echo Setting up compiler environment
call "%IFORT_COMPILER13%\bin\compilervars" intel64
:envexist

erase *.obj erase *.mod
make -j4 -f ..\Makefile intel_win_64_db
