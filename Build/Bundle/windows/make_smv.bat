@echo off

set CURDIR=%CD%

cd ..\..\..\..\smv\Source
git clean -dxf 

cd ..\Build
git clean -dxf

set BUILDDIR=%CD%

call :BUILDLIB LIBS
call :BUILD    background
call :BUILD    dem2fds
call :BUILDFDS fds2ascii
call :BUILD    hashfile
call :BUILD    smokediff
call :BUILD    smokezip
call :BUILD    wind2fds
call :BUILDSMV smokeview
call :BUILD    set_path
call :BUILD    sh2bat
call :BUILD    get_time
echo "bundle complete"

cd %CURDIR%

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
call %script% bot
cd %BUILDDIR%
exit /b /0

:: -------------------------------------------------------------
 :BUILDLIB
:: -------------------------------------------------------------

set dir=%1
set script=make_LIBS_bot

echo
echo ********** building %dir%
echo
cd %CURDIR%\..\..\..\..\Build
cd %dir%\intel_win_64
call %script% bot
cd %BUILDDIR%
exit /b /0

:: -------------------------------------------------------------
 :BUILDSMV
:: -------------------------------------------------------------

set dir=%1
set script=make_%dir%

echo
echo ********** building %dir%
echo
cd %CURDIR%\..\..\..\..\Build
cd %dir%\intel_win_64
call %script% -r bot
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
cd %CURDIR%\..\..\..\..\Build
cd %dir%\intel_win_64
call %script% bot
cd %BUILDDIR%
exit /b /0

:eof
