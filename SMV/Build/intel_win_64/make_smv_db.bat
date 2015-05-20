@echo off
:: setup compiler environment
call ..\..\..\Utilities\Scripts\setup_intel_compilers.bat

set SMV_TESTFLAG=-D pp_BETA
set SMV_TESTSTRING=test_
set KWDIR=..\..\..\Utilities\keyword
set SDIR=..\..\source

erase *.obj *.mod
call %KWDIR%\expand_file %SDIR%\smokeview %SDIR%\shared\string_util.c
make SMV_TESTFLAG="%SMV_TESTFLAG%" SMV_TESTSTRING="%SMV_TESTSTRING%" -f ..\Makefile intel_win_64_db
call %KWDIR%\contract_file %SDIR%\shared\string_util.c
pause
