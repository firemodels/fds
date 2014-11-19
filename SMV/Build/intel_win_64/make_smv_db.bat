@echo off
:: setup compiler environment
call ..\..\..\Utilities\Scripts\setup_intel_compilers.bat

set SMV_TESTFLAG=-D pp_BETA
set SMV_TESTSTRING=test_

erase *.obj *.mod
make SMV_TESTFLAG="%SMV_TESTFLAG%" SMV_TESTSTRING="%SMV_TESTSTRING%" -f ..\Makefile intel_win_64_db
pause
