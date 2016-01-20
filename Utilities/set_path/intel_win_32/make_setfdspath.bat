@echo off
:: setup compiler environment
call ..\..\Scripts\setup_intel_compilers32.bat
if exist "%VS_COMPILER%\vcvars32x86_amd64.bat" call "%VS_COMPILER%\vcvars32x86_amd64"

erase *.obj
erase *.mod
make SHELL="%ComSpec%" -f ..\Makefile intel_win_32
pause
