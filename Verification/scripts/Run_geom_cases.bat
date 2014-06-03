::@echo off

set size=%1

set svn_drive=c:

set SCRIPT_DIR=%CD%
cd %CD%\..
set BASEDIR=%CD%

cd %BASEDIR%\..\
set SVNROOT=%CD%

set TIME_FILE=%SCRIPT_DIR%\smv_case_times.txt

set RUNFDS=call %SVNROOT%\Utilities\Scripts\runfds_win32.bat

:: VVVVVVVVVVVV set parameters VVVVVVVVVVVVVVVVVVVVVV

:: Choose FDS version (size is "", 32 or 64)

if "%size%" == "" (
  set FDSBASE=fds.exe
  set FDSEXE=%FDSBASE%
) else (
  set FDSBASE=fds_%OPENMP%win_%size%.exe
  set FDSEXE=%SVNROOT%\FDS_Compilation\%OPENMP%intel_win_%size%\fds_%OPENMP%win_%size%.exe
)

set BACKGROUNDEXE=%SVNROOT%\Utilities\background\intel_win_32\background.exe

:: Run jobs in background (or not)

set "bg=%BACKGROUNDEXE% -u 85 -d 5 "
:: set bg=

:: ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:: ---------- Ensure that fds exists

call :is_file_installed %FDSEXE%|| exit /b 1

set FDS=%bg%%FDSEXE%
set SH2BAT=%SVNROOT%\Utilities\Data_Processing\sh2bat

echo.
echo FDS=%FDS%
echo.


cd %SCRIPT_DIR%
echo creating FDS case list from SMV_geom_Cases.sh
%SH2BAT% SMV_geom_Cases.sh SMV_geom_Cases.bat

cd %BASEDIR%

echo "smokeview test cases begin" > %TIME_FILE%
date /t >> %TIME_FILE%
time /t >> %TIME_FILE%

call %SCRIPT_DIR%\SMV_geom_Cases.bat

:: erase %SCRIPT_DIR%\SMV_geom_Cases.bat

cd %BASEDIR%
echo "smokeview test cases end" >> %TIME_FILE%
date /t >> %TIME_FILE%
time /t >> %TIME_FILE%

:loop1
tasklist | find /i /c "%FDSBASE%" > temp.out
set /p numexe=<temp.out
echo Number of cases running - %numexe%
if %numexe% == 0 goto finished
Timeout /t 30 >nul 
goto loop1

:finished
echo "FDS geometry cases completed"
goto eof

:: -----------------------------------------
:is_file_installed
:: -----------------------------------------
  
  set program=%1
  %program% -help 1> %SCRIPT_DIR%\exist.txt 2>&1
  type %SCRIPT_DIR%\exist.txt | find /i /c "not recognized" > %SCRIPT_DIR%\count.txt
  set /p nothave=<%SCRIPT_DIR%\count.txt
  if %nothave% GTR 0 (
    echo "***Fatal error: %program% not present"
    echo "Verification suite aborted"
    exit /b 1
  )
  exit /b 0

:eof
