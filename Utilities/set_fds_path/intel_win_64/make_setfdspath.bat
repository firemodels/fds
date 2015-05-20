@echo off
:: setup compiler environment
call ..\..\Scripts\setup_intel_compilers.bat
set KWDIR=..\..\keyword
set SDIR=..\..\..\SMV\source

erase *.obj
call %KWDIR%\expand_file %SDIR%\set_path %SDIR%\shared\string_util.c
make -f ../Makefile intel_win_64
call %KWDIR%\contract_file %SDIR%\shared\string_util.c
pause
