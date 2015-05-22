@echo off
:: setup compiler environment
call ..\..\Scripts\setup_intel_compilers.bat
set KWDIR=..\..\keyword
set SDIR=..\..\..\SMV\source

Title Building getdate for 64 bit Windows

erase *.obj
call %KWDIR%\expand_file %SDIR%\getdate %SDIR%\getdate\main.c
make -f ..\Makefile intel_win_64
call %KWDIR%\contract_file %SDIR%\getdate\main.c
pause

