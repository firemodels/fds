@echo off

Rem windows batch file to build smokeview from the command line

IF "%SETUP_IFORT_COMPILER12%"=="1" GOTO envexist

set SETUP_IFORT_COMPILER12=1

echo Setting up compiler environment
call "%IFORT_COMPILER12%\bin\ifortvars" ia32
call "%IFORT_COMPILER12%\bin\iclvars" ia32
if exist "%VS_COMPILER%\vcvars32x86_amd64.bat" call "%VS_COMPILER%\vcvars32x86_amd64"
:envexist

erase *.obj
erase *.mod
make VPATH="../../../SMV/source/set_path" -f ..\Makefile intel_win_32
pause