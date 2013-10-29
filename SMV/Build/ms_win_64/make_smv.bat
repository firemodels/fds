@echo off

set SMV_TESTFLAG=
set SMV_TESTSTRING=

IF "%COMPILERS_DEFINED%"=="1" GOTO envexist
  set COMPILERS_DEFINED=1
  call "%VS_COMPILER%\..\vcvarsall" x86_amd64
:envexist

if "%1" NEQ "-t" goto endif
  set SMV_TESTFLAG=-D pp_BETA
  set SMV_TESTSTRING=test_
:endif

Rem erase *.obj
make -j4 SMV_TESTFLAG="%SMV_TESTFLAG%" SMV_TESTSTRING="%SMV_TESTSTRING%" -f ..\Makefile ms_win_64
pause
