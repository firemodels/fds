@echo off
:: setup compiler environment
call ..\..\Scripts\setup_intel_compilers32.bat
set KWDIR=..\..\keyword
set SDIR=..\..\..\SMV\source

if exist "%VS_COMPILER%\vcvars32x86_amd64.bat" call "%VS_COMPILER%\vcvars32x86_amd64"
Title Building background for 32 bit Windows

erase *.obj *.mod
call %KWDIR%\expand_file %SDIR%\background %SDIR%\shared\string_util.c
make -f ..\Makefile intel_win_32
call %KWDIR%\contract_file %SDIR%\shared\string_util.c
pause
