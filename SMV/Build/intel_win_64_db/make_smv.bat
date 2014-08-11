@echo off

Rem windows batch file to build smokeview from the command line

IF "%SETUP_IFORT_COMPILER_64%"=="1" GOTO envexist

set SETUP_IFORT_COMPILER_64=1

echo Setting up compiler environment
call "%IFORT_COMPILER14%\bin\compilervars" intel64
:envexist

Title Building debug version of Smokeview for 64 bit Windows

erase *.obj *.mod
make -f ..\Makefile intel_win_64_db
pause

