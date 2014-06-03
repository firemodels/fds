@echo off

set size=%1

echo Creating figures for the Smokeview User's and Verification guides


SETLOCAL
set SCRIPT_DIR=%CD%

cd %CD%\..
set BASEDIR=%CD%

cd %BASEDIR%\..\
set SVNROOT=%CD%

if "%size%" == "" (
  set SMOKEVIEW=smokeview
) else (
  set SMOKEVIEW=%SVNROOT%\SMV\Build\intel_win_%size%\smokeview_win_%size%.exe -bindir %SVNROOT%\SMV\for_bundle
)

set BACKGROUND="background"

call :is_file_installed %BACKGROUND%|| exit /b 1
call :is_file_installed %SMOKEVIEW%|| exit /b 1

set vis="%SVNROOT%\Verification\Visualization"

set RUNFDS=call "%SCRIPT_DIR%\runsmv.bat"
set SH2BAT=%SVNROOT%\Utilities\Data_Processing\sh2bat

echo.
echo converting SMV_geom_Cases.sh case list to SMV_Pictures_Cases.bat

cd %SCRIPT_DIR%
%SH2BAT% SMV_geom_Cases.sh SMV_geom_Pictures_Cases.bat

echo.
echo Generating images

cd %BASEDIR%
call %SCRIPT_DIR%\SMV_geom_Pictures_Cases.bat

cd %SCRIPT_DIR%

erase SMV_Pictures_Cases.bat
erase SMV_DIFF_Pictures_Cases.bat

goto eof

:: -----------------------------------------
  :is_file_installed
:: -----------------------------------------

  set program=%1
  %program% -help 1>> %SCRIPT_DIR%\exist.txt 2>&1
  type %SCRIPT_DIR%\exist.txt | find /i /c "not recognized" > %SCRIPT_DIR%\count.txt
  set /p nothave=<%SCRIPT_DIR%\count.txt
  if %nothave% == 1 (
    echo "***Fatal error: %program% not present"
    echo "Verification aborted"
    exit /b 1
  )
  echo %program% exists
  exit /b 0

:eof
