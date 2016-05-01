@echo off
set release=%1
set from=%2

:: setup compiler environment
if x%from% == xbot goto skip1
call ..\..\..\..\Utilities\Scripts\setup_intel_compilers.bat
:skip1

set SMV_TESTFLAG=
set SMV_TESTSTRING=
if "%release%" == "-r" goto endif
  set SMV_TESTFLAG=-D pp_BETA
  set SMV_TESTSTRING=test_
:endif

erase *.obj *.mod
make -j 4 SHELL="%ComSpec%" SMV_TESTFLAG="%SMV_TESTFLAG%" SMV_TESTSTRING="%SMV_TESTSTRING%" -f ..\Makefile intel_win_64_db > compile.out 2>&1
call :find_smokeview_warnings compile.out

if x%from% == xbot goto skip2
pause
:skip2
goto eof

:: -------------------------------------------------------------
  :find_smokeview_warnings
:: -------------------------------------------------------------

set search_file=%1

grep -v "commands for target" %search_file% > build_warning0.txt
grep -i -A 1 -B 1 remark build_warning0.txt > build_warnings.txt
grep -i -A 1 -B 1 warning build_warning0.txt >> build_warnings.txt
type build_warnings.txt | find /v /c "kdkwokwdokwd"> build_nwarnings.txt
set /p nwarnings=<build_nwarnings.txt
if %nwarnings% GTR 0 (
  echo Warnings:
  echo.
  type build_warnings.txt
)
exit /b

:eof
