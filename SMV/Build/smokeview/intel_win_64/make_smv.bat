@echo off
set release=%1
set from=%2

:: setup compiler environment
if x%from% == xbot goto skip1
call ..\..\..\..\Utilities\Scripts\setup_intel_compilers.bat
:skip1

set SMV_TESTFLAG=
set SMV_TESTSTRING=

Title Building Smokeview for 64 bit Windows
set SMV_TESTFLAG=
set SMV_TESTSTRING=
if "%release%" == "-r" goto endif
  Title Building Test Smokeview for 64 bit Windows
  set SMV_TESTFLAG=-D pp_BETA
  set SMV_TESTSTRING=test_
:endif
set OPT=
if  "x%VS140COMNTOOLS%" == "x" goto endif2
  set OPT=-DHAVE_SNPRINTF -DHAVE_STRUCT_TIMESPEC
:endif2

erase *.obj *.mod
make -j 4 SHELL="%ComSpec%" SMV_TESTFLAG="%SMV_TESTFLAG% %OPT%" SMV_TESTSTRING="%SMV_TESTSTRING%" -f ..\Makefile intel_win_64
if x%from% == xbot goto skip2
pause
:skip2

