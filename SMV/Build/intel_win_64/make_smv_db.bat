@echo off
set arg1=%1

:: setup compiler environment
if x%arg1% == xbot goto skip1
call ..\..\..\Utilities\Scripts\setup_intel_compilers.bat
:skip1

set SMV_TESTFLAG=-D pp_BETA
set SMV_TESTSTRING=test_
set KWDIR=..\..\..\Utilities\keyword
set SDIR=..\..\source

erase *.obj *.mod
make SHELL="%ComSpec%" SMV_TESTFLAG="%SMV_TESTFLAG%" SMV_TESTSTRING="%SMV_TESTSTRING%" -f ..\Makefile intel_win_64_db
if x%arg1% == xbot goto skip2
pause
:skip2
