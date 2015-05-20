@echo off
:: setup compiler environment
call ..\..\Scripts\setup_intel_compilers.bat
set KWDIR=..\..\keyword
set SDIR=..\..\..\SMV\source

Title Building smokediff for 64 bit Windows

erase *.obj *.mod
call %KWDIR%\expand_file %SDIR%\smokediff %SDIR%\shared\string_util.c
make -f ..\Makefile intel_win_64
call %KWDIR%\contract_file %SDIR%\shared\string_util.c
pause

