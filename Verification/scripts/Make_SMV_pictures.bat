@echo off

echo Creating figures for the Smokeview User's and Verification guides


set SCRIPT_DIR=%CD%

cd %CD%\..
set BASEDIR=%CD%

cd %BASEDIR%\..\
set SVNROOT=%CD%

set SMOKEVIEW="smokeview"
set SMOKEZIP="smokezip"
set SMOKEDIFF="smokediff"
set BACKGROUND="background"

set vis="%SVNROOT%\Verification\Visualization"
set wui="%SVNROOT%\Verification\Wui"
set smvug="%SVNROOT%\Manuals\SMV_User_Guide"
set smvvg="%SVNROOT%\Manuals\SMV_Verification_Guide"
set summary="%SVNROOT%\Manuals\Verification_Summary"

set RUNFDS=call "%SCRIPT_DIR%\runsmv.bat"
set RUNCFAST=call "%SCRIPT_DIR%\runsmv.bat"
set SH2BAT=%SVNROOT%\Utilities\Data_Processing\sh2bat

Rem erase summary images

erase %summary%\images\*.png

Rem --------------  user guide ----------------

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

Rem --------------  verification guide ----------------

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

Rem copy images to summary directory

copy %smvug%\SCRIPT_FIGURES\*.png %summary%\images
copy %smvvg%\SCRIPT_FIGURES\*.png %summary%\images

cd %SCRIPT_DIR%
