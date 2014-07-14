"@echo off

:: windows batch file to build smokeview from the command line

IF "%SETUP_IFORT_COMPILER_IA32%"=="1" GOTO envexist

set SETUP_IFORT_COMPILER_IA32=1

echo Setting up compiler environment
call "%IFORT_COMPILER14%\bin\compilervars" ia32
if exist "%VS_COMPILER%\vcvars32x86_amd64.bat" call "%VS_COMPILER%\vcvars32x86_amd64"
:envexist
Title Building smokediff for 32 bit Windows
erase *.obj
erase *.mod
make -f ..\Makefile intel_win_32
pause
