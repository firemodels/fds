@echo off

set release=%1

echo Creating figures for the Smokeview User's and Verification guides


SETLOCAL
set SCRIPT_DIR=%CD%

cd %CD%\..
set BASEDIR=%CD%

cd %BASEDIR%\..\
set SVNROOT=%CD%

if "%release%" == "" (
  set SMOKEDIFF=smokediff
  set SMOKEZIP=smokezip
  set SMOKEVIEW=smokeview
) else (
  set SMOKEDIFF=%SVNROOT%\Utilities\smokediff\intel_win_%release%\smokediff_win_%release%.exe
  set SMOKEVIEW=%SVNROOT%\SMV\Build\intel_win_%release%\smokeview_win_%release%.exe
  set  SMOKEZIP=%SVNROOT%\Utilities\smokezip\intel_win_%release%\smokezip_win_%release%.exe
)

set BACKGROUND="background"

call :is_file_installed %BACKGROUND%|| exit /b 1
call :is_file_installed %SMOKEDIFF%|| exit /b 1
call :is_file_installed %SMOKEVIEW%|| exit /b 1
call :is_file_installed %SMOKEZIP%|| exit /b 1

set vis="%SVNROOT%\Verification\Visualization"
set wui="%SVNROOT%\Verification\Wui"
set smvug="%SVNROOT%\Manuals\SMV_User_Guide"
set smvvg="%SVNROOT%\Manuals\SMV_Verification_Guide"
set summary="%SVNROOT%\Manuals\Verification_Summary"

set RUNFDS=call "%SCRIPT_DIR%\runsmv.bat"
set RUNTFDS=call "%SCRIPT_DIR%\runtsmv.bat"
set RUNWFDS=call "%SCRIPT_DIR%\runsmv.bat"
set RUNCFAST=call "%SCRIPT_DIR%\runsmv.bat"
set SH2BAT=%SVNROOT%\Utilities\Data_Processing\sh2bat

:: erase summary images

erase %summary%\images\*.png

:: --------------  user guide ----------------

cd %smvug%\SCRIPT_FIGURES

echo.
echo erasing Smokeview User guide scripted figures and info files

erase *.png
erase *.help
erase smokeview.version
erase smokediff.version
erase smokezip.version
erase background.version

echo.
echo Creating Smokeview User guide info files

%SMOKEVIEW%  -help > smokeview.help
%SMOKEZIP%   -help > smokezip.help
%SMOKEDIFF%  -help > smokediff.help
%BACKGROUND% -help > background.help

%SMOKEVIEW%  -v > smokeview.version
%SMOKEZIP%   -v > smokezip.version
%SMOKEDIFF%  -v > smokediff.version
%BACKGROUND% -v > background.version

:: --------------  verification guide ----------------

cd %smvvg%\SCRIPT_FIGURES

echo.
echo erasing Smokeview Verification guide scripted figures

erase *.png

echo.
echo converting SMV_Cases.sh case list to SMV_Pictures_Cases.bat

cd %SCRIPT_DIR%
%SH2BAT% SMV_Cases.sh SMV_Pictures_Cases.bat
%SH2BAT% SMV_DIFF_Cases.sh SMV_DIFF_Pictures_Cases.bat

echo.
echo converting plume5c particles to an isosurface

cd %SVNROOT%\Verification\Visualization
%SMOKEZIP% -f -part2iso plumeiso

echo.
echo differencing plume5c and plume5cdelta

%SMOKEDIFF% plume5c plume5cdelta

echo.
echo differencing thouse5 and thouse5delta
%SMOKEDIFF% thouse5 thouse5delta

echo.
echo Generating images

cd %BASEDIR%
call %SCRIPT_DIR%\SMV_Pictures_Cases.bat

cd %BASEDIR%
call %SCRIPT_DIR%\SMV_DIFF_Pictures_Cases.bat

:: copy images to summary directory

echo copying user guide script figures from %smvug%\SCRIPT_FIGURES to %summary%\images
copy %smvug%\SCRIPT_FIGURES\*.png %summary%\images
echo copying verification guide script figures from %smvvg%\SCRIPT_FIGURES to %summary%\images
copy %smvvg%\SCRIPT_FIGURES\*.png %summary%\images
echo copying graysquares figures from %smvvg%\FIGURES to %summary%\images
copy %smvvg%\FIGURES\graysquares.png %summary%\images


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
  exit /b 0

:eof
