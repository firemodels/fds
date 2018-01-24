@echo off

set CURDIR=%CD%

cd ..\..\..\..\smv\Source
git clean -dxf 

cd ..\Build
git clean -dxf

set BUILDDIR=%CD%

call :BUILD LIBS
call :BUILD background
call :BUILD dem2fds
call :BUILDFDS fds2ascii
call :BUILD hashfile
call :BUILD smokediff
call :BUILD wind2fds
call :BUILD smokeview
echo "bundle complete"

goto :eof

:: -------------------------------------------------------------
 :BUILDFDS
:: -------------------------------------------------------------

set dir=%1
set script=make_%dir%

echo
echo ********** building %dir%
echo
cd %CURDIR%\..\..\..\Utilities
cd %dir%\intel_win_64
call %script%
cd %BUILDDIR%
exit /b /0

:: -------------------------------------------------------------
 :BUILD
:: -------------------------------------------------------------

set dir=%1
set script=make_%dir%

echo
echo ********** building %dir%
echo
cd %CURDIR%\..\..\..\Utilities
cd %dir%/intel_win_64
call %script%
cd %BUILDDIR%
exit /b /0

:eof
cd %CURDIR%
