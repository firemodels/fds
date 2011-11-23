@echo off

echo Creating figures for the Smokeview User's and Verification guides


set SCRIPT_DIR=%CD%
set BASEDIR=%CD%\..
set SVNROOT=%BASEDIR%\..\

set SMOKEVIEW="smokeview"
set SMOKEZIP="smokezip_win_64"
set SMOKEDIFF="smokediff_win_64"
set BACKGROUND="background"

set vis="%BASEDIR%\Visualization"
set wui="%BASEDIR%\Wui"
set smvug="%SVNROOT%\Manuals\SMV_User_Guide"
set smvvg="%SVNROOT%\Manuals\SMV_Verification_Guide"
set RUNFDS=call "%SCRIPT_DIR%\runsmv.bat"
set SH2BAT="%SVNROOT%\Utilities\Data_Processing\sh2bat"

cd %smvug%

erase SCRIPT_FIGURES\*.png
erase SCRIPT_FIGURES\*.help
erase SCRIPT_FIGURES\*.version

%SMOKEVIEW% -help > SCRIPT_FIGURES\smokeview.help
%SMOKEVIEW% -version > SCRIPT_FIGURES\smokeview.version
%SMOKEZIP% -help > SCRIPT_FIGURES\smokezip.help
%SMOKEDIFF% -help > SCRIPT_FIGURES\smokediff.help
%SMOKEDIFF% -v > SCRIPT_FIGURES\smokediff.version
%BACKGROUND% -help > SCRIPT_FIGURES\background.help
%BACKGROUND% -version > SCRIPT_FIGURES\background.version

cd %smvvg%
erase SCRIPT_FIGURES\*.version
erase SCRIPT_FIGURES\*.png
%SMOKEVIEW% -version > SCRIPT_FIGURES\smokeview.version

echo converting SMV_Cases.sh case list to SMV_Cases.bat
cd %SCRIPT_DIR%
%SH2BAT% SMV_Cases.sh SMV_Cases.bat

cd %BASEDIR%
call %SCRIPT_DIR%\SMV_Cases.bat

cd %SCRIPT_DIR%
