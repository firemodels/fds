 @echo off
release=%1
set from=%2

:: setup compiler environment
if x%from% == xbot goto skip1
call ..\..\..\Utilities\Scripts\setup_intel_compilers.bat
:skip1

if "%release%" NEQ "-t" goto endif
  set SMV_TESTFLAG=-D pp_BETA
  set SMV_TESTSTRING=test_
:endif

erase *.obj *.mod
make SHELL="%ComSpec%" SMV_TESTFLAG="%SMV_TESTFLAG%" SMV_TESTSTRING="%SMV_TESTSTRING%" -f ..\Makefile intel_win_64_db
if x%from% == xbot goto skip2
pause
:skip2
