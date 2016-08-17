@echo off
:: setup MS compiler environment
::call "%VS_COMPILER%\vcvarsall.bat" amd64

:: setup intel compiler environment
call ..\..\..\Utilities\Scripts\setup_intel_compilers.bat


set SMV_TESTFLAG=
set SMV_TESTSTRING=

Title Building Smokeview for 64 bit Windows
if "%1" NEQ "-t" goto endif
  Title Building Test Smokeview for 64 bit Windows
  set SMV_TESTFLAG=-D pp_BETA
  set SMV_TESTSTRING=test_
:endif

:: erase *.obj *.mod
make SMV_TESTFLAG="%SMV_TESTFLAG%" SMV_TESTSTRING="%SMV_TESTSTRING%" -f ..\Makefile ms_win_64
pause

