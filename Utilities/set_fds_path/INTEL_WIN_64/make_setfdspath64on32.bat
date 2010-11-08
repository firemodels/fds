@echo off

Rem windows batch file to build setfdspath

IF "%SETUP_IFORT_COMPILER11%"=="1" GOTO envexist

set SETUP_IFORT_COMPILER11=1

echo Setting up compiler environment
call "%IFORT_COMPILER11%\bin\ifortvars" ia32_intel64
call "%IFORT_COMPILER11%\bin\iclvars" ia32_intel64

:envexist
erase *.obj
make VPATH="../../../SMV_5/source/set_path" -f ../Makefile intel_win_64
pause